using System.Collections;
using System.Collections.Generic;
using Unity.Burst;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Jobs;
using UnityEngine;
using UnityEngine.Jobs;
using UnityEngine.Rendering;
using GluonGui.WorkspaceWindow.Views.WorkspaceExplorer.Explorer;
using UnityEngine.UIElements;



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
        [SerializeField] private float minComputeTimeChunks = 0.025f;

        [Header("Physics")]
        [SerializeField] private bool usePrimitiveColliders = false;
        [SerializeField] private float primitiveSizeFactor = 1.0f;
        [SerializeField] private bool primitivePreferHoles = true;
        [SerializeField] private PhysicMaterial phyMat = null;
        [SerializeField] private float voxelMass = 0.1f;
        [SerializeField] private RigidbodyInterpolation interpolation = RigidbodyInterpolation.Interpolate;

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
            if (IsRuntimeCreated() == true) return;
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

        private bool IsRuntimeCreated()
        {

            return Application.isPlaying == true && gameObject.GetInstanceID() < 0;
        }

        private void OnValidateDelayed()
        {
            if (Application.isPlaying == true)
            {
                if (IsRuntimeCreated() == true) return;
                Debug.Log(transform.name + " regenerated MeshGen by OnValidate");
            }

            Dispose();
            Init();
        }
#endif

        private bool isInitilized = false;
        private MeshCollider col = null;
        private MeshFilter mf = null;
        private MeshRenderer mr = null;
        private Mesh mesh = null;
        private Rigidbody rb = null;
        private bool rbWasKin = false;

        private void SetMO(MarchingObject newMO)
        {
            cm_job = new()
            {
                mo = newMO,
                result = new(Allocator.Persistent)
            };

            meshNeedsComputing = true;
        }

        private void Init()
        {
            if (isInitilized == true || tempDisableInit == true) return;
            isInitilized = true;
            lastShapeID = _currentShapeID;
            lastConfigID = _currentConfigID;

            if (cm_job.result.IsCreated == false)
            {
                rb = GetComponent<Rigidbody>();
                col = this.GetOrAddComponent<MeshCollider>(out _);
                col.sharedMesh = startShape;
                col.enabled = true;

                if (startShape != null) NewMarchingObject(col);
                else if (TryGetComponent(out Collider othCol) == true && othCol.enabled == true && othCol != col)
                {
                    NewMarchingObject(othCol);
                }
                else
                {
                    isInitilized = false;
                    return;
                }

                void NewMarchingObject(Collider col)
                {
                    SetMO(new(col, kinematicVoxelsOverlap, tryFillConvaveInteriors, tryStripOuterLayer, surfaceLevel, isoLevel));
                }
            }

            mesh = new();
            mesh.MarkDynamic();
            if (col != null)
            {
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
            }

            if (rb != null)
            {
                rb.interpolation = interpolation;
                rbWasKin = rb.isKinematic;
                UpdateRigidbody();
            }

            mf = this.GetOrAddComponent<MeshFilter>(out _);
            mr = this.GetOrAddComponent<MeshRenderer>(out _);
            mf.sharedMesh = mesh;
            mr.sharedMaterial = mat;

            chunkBaseName = transform.name + "Chunk_";
            chunkBaseLayer = gameObject.layer;
            primBaseLayer = gameObject.layer;
            primBaseName = transform.name + "Prim_";
            primBaseRadius = _mo.voxSize.Average() / (2 / primitiveSizeFactor);
            pendingActionsA = new(8, Allocator.Persistent);
            pendingActionsB = new(8, Allocator.Persistent);

            bc_job = new()
            {
                meshID = mesh.GetInstanceID(),
                convex = col != null && col.convex,
            };

            allMeshGenObjects.Add(this);
            meshNeedsComputing = true;
            timeSinceStartedComputing = 69420.0f;
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            ComputeMesh_end();
            BakeCollision_end();
            allMeshGenObjects.Remove(this);
            DisposePrims();
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
        private int ignoreFrames = 0;

        private void Update()
        {
            if (isInitilized == false) return;
            if (usePrimitiveColliders == true) UpdatePrims();
            if (timeSinceStartedComputing < minComputeTime)
            {
                timeSinceStartedComputing += Time.deltaTime;
                return;
            }

            if (ignoreFrames > 0)
            {
                ignoreFrames--;
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

            timeSinceStartedComputing = chunksOnly == false ? 0.0f : (minComputeTime - minComputeTimeChunks);
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
                    if (usePrimitiveColliders == true)
                    {
                        foreach (int voxI in _mo.addedRemovedVoxsI)
                        {
                            if (voxI < 0)
                            {//Remove colliders instantly
                                int vI = Mathf.Abs(voxI);
                                voxsToUpdate.Remove(vI);//If added and deleted before processed
                                if (voxIToCol.Remove(vI, out SphereCollider col) == false) continue;
                                col.gameObject.SetActive(false);
                                inactiveCols.Enqueue(col);
                                continue;
                            }

                            voxsToUpdate.Add(voxI);
                        }
                    }

                    UpdateRigidbody();
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

        private static readonly Queue<SphereCollider> inactiveCols = new(64);
        private static readonly HashSet<MeshGenObject> allMeshGenObjects = new(16);
        private readonly Dictionary<int, SphereCollider> voxIToCol = new(128);
        private string primBaseName = "meshGenObj_";
        private int primBaseLayer = 0;
        private float primBaseRadius = MeshGenGlobals.voxelSizeWorld;

        private void DisposePrims()
        {
            bool playing = Application.isPlaying;
            foreach (SphereCollider col in voxIToCol.Values)
            {
                if (col == null) continue;
                if (playing == false) DestroyImmediate(col.gameObject);
                else Destroy(col.gameObject);
            }

            voxIToCol.Clear();
            if (allMeshGenObjects.Count > 0) return;

            foreach (SphereCollider col in inactiveCols)
            {
                if (col == null) continue;
                if (playing == false) DestroyImmediate(col.gameObject);
                else Destroy(col.gameObject);
            }

            inactiveCols.Clear();
            return;
        }

        private readonly HashSet<int> voxsToUpdate = new(128);
        private readonly int[] voxsUpdated = new int[maxPrimUpdatesPerFrame];
        private const int maxPrimUpdatesPerFrame = 20;//Creating 20 takes ~0.9ms in Editor

        private unsafe void UpdatePrims()
        {
            if (voxsToUpdate.Count == 0) return;
            int loopCount = -1;

            foreach (int voxI in voxsToUpdate)
            {
                if (++loopCount >= maxPrimUpdatesPerFrame) break;
                voxsUpdated[loopCount] = voxI;

                if (voxIToCol.TryGetValue(voxI, out SphereCollider col) == false)
                {
                    bool foundCol = false;

                    while (inactiveCols.TryDequeue(out col) == true)
                    {
                        if (col == null) continue;
                        if (col.transform.parent != transform)//Worth checking?
                        {
                            col.transform.SetParent(transform, false);
                            col.transform.localScale = Vector3.one;
                            col.gameObject.layer = primBaseLayer;
                            col.sharedMaterial = phyMat;
                            col.radius = primBaseRadius;
                        }

                        col.gameObject.SetActive(true);
                        foundCol = true;
                        break;
                    }

                    if (foundCol == false)
                    {
                        GameObject newO = new(primBaseName + voxI, typeof(SphereCollider));
                        newO.transform.SetParent(transform, false);
                        newO.transform.localScale = Vector3.one;
                        newO.layer = primBaseLayer;
                        newO.hideFlags = HideFlags.DontSave;
                        col = newO.GetComponent<SphereCollider>();
                        col.sharedMaterial = phyMat;
                        col.radius = primBaseRadius;
                    }

                    voxIToCol.Add(voxI, col);
                }
                else if (primitivePreferHoles == false) continue;//Only voxState changed

                if (primitivePreferHoles == true)
                {
                    int vI = voxI;
                    byte vState = _mo.voxsState[voxI];//Can read while burst is writing, potential race, who cares
                    if ((vState & (1 << 1)) != 0) vI += _mo.vCountYZ;
                    if ((vState & (1 << 2)) != 0) vI += _mo.vCountZ;
                    if ((vState & (1 << 3)) != 0) vI += 1;

                    col.transform.localPosition = _mo.VoxIndexToPosL(vI);
                    continue;
                }

                col.transform.localPosition = _mo.VoxIndexToPosL(voxI)
                    + (_mo.voxSize * 0.5f);
            }

            if (loopCount >= maxPrimUpdatesPerFrame)
            {
                foreach (int voxI in voxsUpdated)
                {
                    voxsToUpdate.Remove(voxI);
                }
            }
            else voxsToUpdate.Clear();
        }

        private void UpdateRigidbody()
        {
            if (rb == null) return;
            int solidCount = _mo.solidVoxCount.Value;
            if (solidCount <= 0)
            {
                rb.isKinematic = true;
                return;
            }

            rb.isKinematic = rbWasKin;
            rb.mass = voxelMass * solidCount;
        }

        #endregion Collision

        #region Chunks

        private string chunkBaseName = "chunk_";
        private int chunkBaseLayer = 0;
        private int nextChunkCount = 0;
        private static bool tempDisableInit = false;

        public unsafe void TryCreateChunkFrom(ref MarchingObject.Unmanaged mo)
        {//~0.25ms total avg (Editor) 
            if (mo.overlappingVoxCount <= 0) return;
            meshNeedsComputing = true;
            tempDisableInit = true;

            //Maintain voxCount, needed for rbMass
            GameObject newO = new(chunkBaseName + nextChunkCount++, typeof(MeshGenObject),
                typeof(Rigidbody), typeof(MeshFilter), typeof(MeshRenderer));
            newO.transform.SetParent(transform.parent);
            newO.transform.localScale = transform.localScale;
            transform.GetLocalPositionAndRotation(out Vector3 pos, out Quaternion rot);
            newO.transform.SetLocalPositionAndRotation(pos, rot);
            newO.layer = chunkBaseLayer;

            Rigidbody newRb = newO.GetComponent<Rigidbody>();
            newRb.interpolation = interpolation;
            if (this.rb != null)
            {
                newRb.velocity = this.rb.velocity;
                newRb.angularVelocity = this.rb.angularVelocity;
                newRb.useGravity = this.rb.useGravity;
            }

            MeshGenObject mgo = newO.GetComponent<MeshGenObject>();
            mgo.SetMO(new(mo));
            mgo.usePrimitiveColliders = true;
            mgo.interpolation = interpolation;
            mgo.mat = mat;
            mgo.phyMat = phyMat;
            mgo.isoLevel = isoLevel;
            mgo.surfaceLevel = surfaceLevel;
            mgo.minComputeTime = minComputeTime;
            mgo.primitivePreferHoles = primitivePreferHoles;
            mgo.primitiveSizeFactor = primitiveSizeFactor;
            mgo.voxelMass = voxelMass;
            mgo.rb = newRb;

            tempDisableInit = false;
            mgo.Init();//~30% of creation time is this ~0.085ms (Editor)

            ignoreFrames = 2;
            mo.overlappingVoxCount = -1;
        }

        #endregion Chunks
    }
}


