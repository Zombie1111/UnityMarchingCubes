using System.Collections;
using System.Collections.Generic;
using Unity.Burst;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Jobs;
using UnityEngine;
using UnityEngine.Assertions.Must;
using UnityEngine.Rendering;

namespace zombGen
{
    [ExecuteAlways]
    public class MeshGenObject : MonoBehaviour
    {
        [Header("Setup")]
        [Tooltip("If null, tries to use collider attatch to this trans as shape")]
        [SerializeField] private Mesh startShape = null;
        [SerializeField] private bool tryFillConvaveInteriors = true;
        [SerializeField] private bool tryStripOuterLayer = false;
        [SerializeField] private float surfaceLevel = 0.0f;
        [SerializeField] private float isoLevel = 1.0f;

        [Header("Runtime")]
        [SerializeField] private float minComputeTime = 0.05f;

        [Header("Debug")]
        [SerializeField] private bool drawStartShape = false;

#if UNITY_EDITOR
        private void OnDrawGizmosSelected()
        {
            if (drawStartShape == true && startShape != null)
            {
                Gizmos.DrawMesh(startShape, 0, transform.position, transform.rotation, transform.lossyScale);
            }
        }
#endif

        private void OnEnable()
        {
            Init();
        }

        private void OnDisable()
        {
            Dispose();
        }

        private int lastShapeID = 0;
        private int _currentShapeID => startShape != null ? startShape.GetInstanceID() : 0;
        private int lastConfigID = 0;
        private int _currentConfigID
        {
            get
            {
                int id = 0;
                if (tryFillConvaveInteriors == true) id += 9823;
                if (tryStripOuterLayer == true) id += 65;
                id += (int)(surfaceLevel * 1000.0f);
                id += (int)(isoLevel * 1000.0f);
                return id;
            }
        }


#if UNITY_EDITOR
        private void OnValidate()
        {
            if (lastShapeID == _currentShapeID
                && lastConfigID == _currentConfigID) return;

            if (Application.isPlaying == true)
                Debug.Log(transform.name + " regenerated MeshGen by OnValidate");

            Dispose();
            Init();
        }
#endif

        private bool isInitilized = false;
        private MeshCollider col = null;
        private MeshFilter mf = null;
        private MeshRenderer mr = null;
        private Mesh mesh = null;

        private void Init()
        {
            if (isInitilized == true) return;
            isInitilized = true;
            lastShapeID = _currentShapeID;
            lastConfigID = _currentConfigID;
            col = this.GetOrAddComponent<MeshCollider>(out _);
            mf = this.GetOrAddComponent<MeshFilter>(out _);
            mr = this.GetOrAddComponent<MeshRenderer>(out _);
            col.sharedMesh = startShape;
            MarchingObject mo;

            if (startShape != null) mo = NewMarchingObject(col);
            else if (TryGetComponent(out Collider othCol) == true && othCol.enabled == true)
            {
                mo = NewMarchingObject(othCol);
            }
            else
            {
                isInitilized = false;
                return;
            }

            MarchingObject NewMarchingObject(Collider col)
            {
                return new(col, tryFillConvaveInteriors, tryStripOuterLayer, surfaceLevel, isoLevel);
            }

            mesh = new();
            mesh.MarkDynamic();
            col.sharedMesh = mesh;
            mf.sharedMesh = mesh;

            pendingActionsA = new(8, Allocator.Persistent);
            pendingActionsB = new(8, Allocator.Persistent);

            cm_job = new()
            {
                mo = mo,
                result = new(Allocator.Persistent)
            };

            bm_job = new()
            {
                meshID = mesh.GetInstanceID(),
                convex = col.convex,
            };

            meshNeedsComputing = true;
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            isInitilized = false;
            ComputeMesh_end();
            BakeMesh_end();
            cm_job.mo.Dispose();
            cm_job.result.Dispose();
            pendingActionsA.Dispose();
            pendingActionsB.Dispose();
            if (col != null) col.sharedMesh = null;
            mesh.Clear();
        }
        
        private const MeshUpdateFlags _defaultUpdateFlags =
              MeshUpdateFlags.DontRecalculateBounds
            | MeshUpdateFlags.DontValidateIndices
            | MeshUpdateFlags.DontNotifyMeshUsers
            | MeshUpdateFlags.DontResetBoneBounds;

        private const MeshColliderCookingOptions _defaultCookingFlags =
            MeshColliderCookingOptions.CookForFasterSimulation
            | MeshColliderCookingOptions.EnableMeshCleaning
            | MeshColliderCookingOptions.WeldColocatedVertices
            | MeshColliderCookingOptions.UseFastMidphase;

        private bool meshNeedsComputing = true;

        private void Update()
        {
            if (isInitilized == false) return;
            if (timeSinceStartedComputing < minComputeTime)
            {
                timeSinceStartedComputing += Time.deltaTime;
                return;
            }

            ComputeMesh_end();
            if (meshNeedsComputing == true) ComputeMesh_start();
        }

        private void LateUpdate()
        {
            if (isInitilized == false) return;
            BakeMesh_end();
        }

        private bool flipped = false;
        private NativeList<Action> pendingActionsA;
        private NativeList<Action> _pendingActionsWrite => flipped == false ? pendingActionsA : pendingActionsB;
        private NativeList<Action> pendingActionsB;
        private NativeList<Action> _pendingActionsRead => flipped == false ? pendingActionsB : pendingActionsA;

        private void Flip()
        {
            flipped = !flipped;
            _pendingActionsWrite.Clear();
        }

        public void PerformAction(Vector3 pos, float radius, Action.Type type)
        {
            _pendingActionsWrite.Add(new(pos, radius, type));
            meshNeedsComputing = true;
        }

        #region ComputeMesh

        private ComputeMesh cm_job;
        private JobHandle cm_handle;
        private bool cm_isActive = false;
        private float timeSinceStartedComputing = 69420.0f;

        private void ComputeMesh_start()
        {
            if (cm_isActive == true) return;
            cm_isActive = true;

            Flip();
            cm_job.wToL = transform.worldToLocalMatrix;
            cm_job.mda = Mesh.AllocateWritableMeshData(1);
            cm_job.actions = _pendingActionsRead.AsArray().AsReadOnly();
            cm_handle = bm_isActive == true ? cm_job.Schedule(bm_handle) : cm_job.Schedule();
            timeSinceStartedComputing = 0.0f;
            meshNeedsComputing = false;
        }

        private void ComputeMesh_end()
        {
            if (cm_isActive == false) return;
            cm_handle.Complete();
            cm_isActive = false;

            ComputeMesh.Result result = cm_job.result.Value;
            Mesh.ApplyAndDisposeWritableMeshData(cm_job.mda, mesh, _defaultUpdateFlags);
            mesh.bounds = result.meshBoundsLocal;

            BakeMesh_start();
        }

        private BakeMesh bm_job;
        private JobHandle bm_handle;
        private bool bm_isActive = false;

        public readonly struct Action
        {
            public readonly Vector3 posW;
            public readonly float radiusW;
            public readonly Type type;

            public Action(Vector3 pos, float radius, Type type = Type.remove)
            {
                posW = pos;
                radiusW = radius;
                this.type = type;
            }

            public enum Type
            {
                remove,
                create
            }
        }

        [BurstCompile]
        private struct ComputeMesh : IJob
        {
            public readonly struct Result
            {
                public Result(Bounds boundsL)
                {
                    meshBoundsLocal = boundsL;
                }

                public readonly Bounds meshBoundsLocal;
            }

            [NativeDisableUnsafePtrRestriction] public MarchingObject mo;
            public Mesh.MeshDataArray mda;
            public Matrix4x4 wToL;
            public NativeReference<Result> result;
            public NativeArray<Action>.ReadOnly actions;

            public void Execute()
            {
                float avgScale = wToL.lossyScale.Average();

                foreach (Action a in actions)
                {
                    int voxI = mo.PosLToVoxIndex(wToL.MultiplyPoint3x4(a.posW));
                    mo.SetVoxelsAround(voxI, a.radiusW, a.type == Action.Type.remove ? -1 : 1);
                }

                mo.Meshify(mda[0], out Bounds bl);
                result.Value = new(bl);
            }
        }

        private void BakeMesh_start()
        {
            if (bm_isActive == true || isInitilized == false) return;

            bm_isActive = true;
            bm_handle = bm_job.Schedule();
        }

        private void BakeMesh_end()
        {
            if (bm_isActive == false) return;
            bm_isActive = false;

            bm_handle.Complete();
            //mesh.MarkModified();//Does not do a shit
            col.sharedMesh = null;
            col.sharedMesh = mesh;//Seems to be the only way to get the collider to actually update
        }

        [BurstCompile]
        private struct BakeMesh : IJob
        {
            public int meshID;
            public bool convex;

            public void Execute()
            {
                Physics.BakeMesh(meshID, convex, _defaultCookingFlags);
            }
        }

        #endregion ComputeMesh
    }
}


