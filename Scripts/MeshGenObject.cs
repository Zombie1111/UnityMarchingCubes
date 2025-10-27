using System.Collections;
using System.Collections.Generic;
using Unity.Burst;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Jobs;
using UnityEngine;
using UnityEngine.Jobs;
using UnityEngine.Rendering;
using Codice.Client.BaseCommands.Annotate;


#if UNITY_EDITOR
using UnityEditor;
#endif

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
        [Tooltip("Voxels overlapping with this will become kinematic")]
        [SerializeField] private Collider kinematicVoxelsOverlap = null;

        [Header("Generic")]
        [SerializeField] private Material mat = null;
        [SerializeField] private float minComputeTime = 0.05f;

        [Header("Physics")]
        [SerializeField] private bool usePrimitiveColliders = false;
        [SerializeField] private float primitiveSizeFactor = 1.0f;
        [SerializeField] private bool primitivePreferHoles = true;
        [SerializeField] private PhysicMaterial phyMat = null;
        [SerializeField] private float voxelMass = 0.1f;

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
                if (usePrimitiveColliders == true) id += 545;
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

            EditorApplication.delayCall += OnValidateDelayed;
        }

        private void OnValidateDelayed()
        {
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
            col.sharedMesh = startShape;
            col.enabled = true;

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
                return new(col, kinematicVoxelsOverlap, tryFillConvaveInteriors, tryStripOuterLayer, surfaceLevel, isoLevel);
            }

            mesh = new();
            mesh.MarkDynamic();
            if (usePrimitiveColliders == false)
            {
                col.sharedMesh = mesh;
                col.sharedMaterial = phyMat;
            }
            else
            {
                col.sharedMesh = null;
                col.enabled = false;
            }

            mf = this.GetOrAddComponent<MeshFilter>(out _);
            mr = this.GetOrAddComponent<MeshRenderer>(out _);
            mf.sharedMesh = mesh;
            mr.sharedMaterial = mat;

            chunkBaseName = transform.name + "Chunk_";
            primBaseLayer = gameObject.layer;
            primBaseName = transform.name + "Prim_";
            primBaseRadius = mo.voxSize.Average() / (2 / primitiveSizeFactor);
            pendingActionsA = new(8, Allocator.Persistent);
            pendingActionsB = new(8, Allocator.Persistent);

            cm_job = new()
            {
                mo = mo,
                result = new(Allocator.Persistent)
            };

            bc_job = new()
            {
                meshID = mesh.GetInstanceID(),
                convex = col.convex,
            };

            meshNeedsComputing = true;
            if (usePrimitiveColliders == true)
                OnVoxelsChanged += DoOnVoxelsChanged;
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            ComputeMesh_end();
            BakeCollision_end();
            OnVoxelsChanged?.Invoke(ChangeType.disposing, _mo.addedRemovedVoxsI);
            OnVoxelsChanged -= DoOnVoxelsChanged;
            isInitilized = false;

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
        private bool chunksNeedsComputing = true;

        private void Update()
        {
            if (isInitilized == false) return;
            if (timeSinceStartedComputing < minComputeTime)
            {
                timeSinceStartedComputing += Time.deltaTime;
                return;
            }

            ComputeMesh_end();
            if (meshNeedsComputing == true || chunksNeedsComputing == true) ComputeMesh_start();
        }

        private void LateUpdate()
        {
            if (isInitilized == false) return;
            BakeCollision_end();
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

        public enum ChangeType
        {
            removedCreated,
            disposing,
        }

        public delegate void Event_OnVoxelsChanged(ChangeType type, NativeList<int> voxsRemovedCreated);
        /// <summary>
        /// if removedCreated: voxsRemovedCreated contains voxel indexs, negative if removed (Only read durring callback)
        /// if initilized: all voxels currently solid was created
        /// if disposing: all voxels currently solid will be removed
        /// </summary>
        public event Event_OnVoxelsChanged OnVoxelsChanged;

        public MarchingObject _mo => cm_job.mo;

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
            bool chunksOnly = chunksNeedsComputing == true && meshNeedsComputing == false;

            if (chunksOnly == false)
            {
                Flip();
                cm_job.wToL = transform.worldToLocalMatrix;
                cm_job.mda = Mesh.AllocateWritableMeshData(1);
                cm_job.actions = _pendingActionsRead.AsArray().AsReadOnly();
            }
            else cm_job.mda = Mesh.AllocateWritableMeshData(0);

            cm_job.chunksOnly = chunksOnly;
            cm_handle = bc_isActive == true ? cm_job.Schedule(bc_handle) : cm_job.Schedule();

            timeSinceStartedComputing = chunksOnly == false ? 0.0f : (minComputeTime / 2);
            meshNeedsComputing = false;
            chunksNeedsComputing = false;
        }

        private void ComputeMesh_end()
        {
            if (cm_isActive == false) return;
            cm_handle.Complete();
            cm_isActive = false;
            bool wasChunksOnly = cm_job.chunksOnly;

            if (wasChunksOnly == false)
            {
                ComputeMesh.Result result = cm_job.result.Value;
                Mesh.ApplyAndDisposeWritableMeshData(cm_job.mda, mesh, _defaultUpdateFlags);
                mesh.bounds = result.meshBoundsLocal;

                if (_mo.addedRemovedVoxsI.Length > 0)
                {
                    OnVoxelsChanged?.Invoke(ChangeType.removedCreated, _mo.addedRemovedVoxsI);
                    _mo.addedRemovedVoxsI.Clear();
                }
            }
            else cm_job.mda.Dispose();
            
            if (_mo.voxsToCheck.Length > 0)
                chunksNeedsComputing = true;

            MarchingObject.Unmanaged mo = _mo.newMO.Value;
            TryCreateChunkFrom(ref mo);
            cm_job.mo.newMO.Value = mo;
            if (wasChunksOnly == false) BakeCollision_start();
        }

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
            public bool chunksOnly;

            public void Execute()
            {
                if (chunksOnly == false)
                {
                    foreach (Action a in actions)
                    {
                        int voxI = mo.PosLToVoxIndex(wToL.MultiplyPoint3x4(a.posW));
                        mo.SetVoxelsAround(voxI, a.radiusW, a.type == Action.Type.remove ? -1 : 1);
                    }

                    mo.Meshify(mda[0], out Bounds bl);
                    result.Value = new(bl);
                }

                mo.TryGetChunk();
            }
        }

        #endregion ComputeMesh

        #region Collision

        private BakeCollision bc_job;
        private JobHandle bc_handle;
        private bool bc_isActive = false;

        private void BakeCollision_start()
        {
            if (bc_isActive == true || isInitilized == false || usePrimitiveColliders == true) return;

            bc_isActive = true;
            bc_handle = bc_job.Schedule();
        }

        private void BakeCollision_end()
        {
            if (bc_isActive == false) return;
            bc_isActive = false;

            bc_handle.Complete();
            //mesh.MarkModified();//Does not do a shit
            col.sharedMesh = null;
            col.sharedMesh = mesh;//Seems to be the only way to get the collider to actually update
        }

        [BurstCompile]
        private struct BakeCollision : IJob
        {
            public int meshID;
            public bool convex;

            public void Execute()
            {
                Physics.BakeMesh(meshID, convex, _defaultCookingFlags);
            }
        }

        private readonly Queue<SphereCollider> inactiveCols = new(16);
        private readonly Dictionary<int, SphereCollider> voxIToCol = new(64);
        private string primBaseName = "meshGenObj_";
        private int primBaseLayer = 0;
        private float primBaseRadius = MeshGenGlobals.voxelSizeWorld;

        private unsafe void DoOnVoxelsChanged(ChangeType type, NativeList<int> voxsRemovedCreated)
        {
            if (type == ChangeType.disposing)
            {
                bool playing = Application.isPlaying;
                foreach (SphereCollider col in inactiveCols)
                {
                    if (playing == false) DestroyImmediate(col.gameObject);
                    else Destroy(col.gameObject);
                }

                foreach (SphereCollider col in voxIToCol.Values)
                {
                    if (col == null) continue;
                    if (playing == false) DestroyImmediate(col.gameObject);
                    else Destroy(col.gameObject);
                }

                voxIToCol.Clear();
                inactiveCols.Clear();
                return;
            }

            int vCount = voxsRemovedCreated.Length;
            for (int i = 0; i < vCount; i++)
            {
                int vI = voxsRemovedCreated[i];
                bool created = vI >= 0;
                vI = Mathf.Abs(vI);
                SphereCollider col;

                if (created == true)
                {
                    if (voxIToCol.TryGetValue(vI, out col) == false)
                    {
                        if (inactiveCols.TryDequeue(out col) == false)
                        {
                            GameObject newO = new(primBaseName + vI);
                            newO.transform.SetParent(transform, false);
                            newO.transform.localScale = Vector3.one;
                            newO.layer = primBaseLayer;
                            newO.hideFlags = HideFlags.DontSave;
                            col = newO.AddComponent<SphereCollider>();
                            col.sharedMaterial = phyMat;
                            col.radius = primBaseRadius;
                        }
                        else
                        {
                            col.gameObject.SetActive(true);
                        }

                        voxIToCol.Add(vI, col);
                    }
                    else if (primitivePreferHoles == false) continue;//Only voxState changed

                    if (primitivePreferHoles == true)
                    {
                        int vI2 = vI;
                        byte vState = _mo.voxsState[vI];
                        if ((vState & (1 << 1)) != 0) vI2 += _mo.vCountYZ;
                        if ((vState & (1 << 2)) != 0) vI2 += _mo.vCountZ;
                        if ((vState & (1 << 3)) != 0) vI2 += 1;

                        col.transform.localPosition = _mo.VoxIndexToPosL(vI2);
                    }
                    else
                    {
                        col.transform.localPosition = _mo.VoxIndexToPosL(vI)
                            + (_mo.voxSize * 0.5f);
                    }

                    continue;
                }

                if (voxIToCol.Remove(vI, out col) == false) continue;
                col.gameObject.SetActive(false);
                inactiveCols.Enqueue(col);
            }
        }

        #endregion Collision

        #region Chunks

        private string chunkBaseName = "chunk_";
        private int nextChunkCount = 0;

        public unsafe void TryCreateChunkFrom(ref MarchingObject.Unmanaged mo)
        {
            if (mo.overlappingVoxCount <= 0) return;

            Debug.Log("Got new chunk " + mo.vCountXYZ + " " + _mo.vCountXYZ);
            meshNeedsComputing = true;

            //Maintain voxCount, needed for rbMass
            //GameObject newO = new(chunkBaseName + nextChunkCount++);
            //Rigidbody rb = newO.AddComponent<Rigidbody>();
            //rb.mass = voxelMass * mo.overlappingVoxCount;
            //
            //MeshGenObject mgo = newO.AddComponent<MeshGenObject>();

            mo.overlappingVoxCount = -1;
        }

        #endregion Chunks
    }
}


