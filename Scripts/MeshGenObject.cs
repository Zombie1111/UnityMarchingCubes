//Copyright 2025 David Westberg (MIT) https://github.com/Zombie1111/UnityMarchingCubes
using System.Collections.Generic;
using Unity.Burst;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Jobs;
using UnityEngine;
using UnityEngine.Jobs;
using UnityEngine.Rendering;
using Unity.Mathematics;

#if UNITY_EDITOR
using UnityEditor;
#endif

namespace zombGen
{
    [ExecuteAlways]
    public class MeshGenObject : MonoBehaviour
    {
        [Header("Generation")]
        [Tooltip("If null, tries to use collider attatch to this trans as shape")]
        [SerializeField] private Mesh startShape = null;
        [SerializeField] private bool tryFillConvaveInteriors = true;
        [SerializeField] private bool tryStripOuterLayer = false;
        [SerializeField] private float surfaceLevel = 0.3f;
        [SerializeField] private float isoLevel = 0.7f;
        [Tooltip("Voxels overlapping with this will become kinematic")]
        [SerializeField] private Collider kinematicVoxelsOverlap = null;
        [SerializeField] private bool smoothNormals = true;

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

        [Header("Breaking")]
        [SerializeField] private BreakMode breakMode = BreakMode.rigidbody;
        [SerializeField] private float forcePerBreakRadius = 40.0f;
        [SerializeField] private int maxBreakRadius = 5;
        [SerializeField] private float voxRemoveInterval = 0.01f;
        [SerializeField] private float perVoxIntervalFactor = 0.03f;
        //[SerializeField] private AudioReference breakAudio = new();
        //private readonly AudioProps breakAudioProps = new();

        public enum BreakMode
        {
            none,
            rigidbody,
            remove,
        }

        [Header("Debug")]
        [SerializeField] private bool drawStartShape = false;

        #region Main
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
                id += (int)breakMode * 46;
                return id;
            }
        }

#if UNITY_EDITOR
        private void OnValidate()
        {
            if (isChunk == true || tempDisableInit == true) return;
            EditorApplication.delayCall += OnValidateDelayed;
        }

        private bool isChunk = false;

        private void OnValidateDelayed()
        {
            if (lastShapeID == _currentShapeID
            && lastConfigID == _currentConfigID)
            {
                ComputeMesh_end();
                if (col != null) col.sharedMaterial = phyMat;
                if (mr != null) mr.sharedMaterial = mat;
                cm_job.smoothNormals = smoothNormals;
                meshNeedsComputing = true;
                return;
            }


            if (Application.isPlaying == true)
            {
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
                result = new(Allocator.Persistent),
                smoothNormals = smoothNormals,
                breakMode = breakMode,
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
                else if (TryGetComponent(out Collider othCol) == true && othCol != col)
                {
                    bool wasEnabled = othCol.enabled;
                    bool wasTrigger = othCol.isTrigger;
                    othCol.isTrigger = false;
                    othCol.enabled = true;
                    NewMarchingObject(othCol);
                    othCol.enabled = wasEnabled;
                    othCol.isTrigger = wasTrigger;
                }
                else
                {
                    isInitilized = false;
                    return;
                }

                void NewMarchingObject(Collider col)
                {
                    SetMO(new(col, kinematicVoxelsOverlap, tryFillConvaveInteriors, tryStripOuterLayer,
                        surfaceLevel, isoLevel, usePrimitiveColliders, breakMode));
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
            primBaseRadius = _mo.data.voxSize.Average() / (2 / primitiveSizeFactor);
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
            //InitVoxels();
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            ComputeMesh_end();
            BakeCollision_end();
            allMeshGenObjects.Remove(this);
            DisposePrims();
            //DisposeVoxels();
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
            BakeCollision_end();
            if (usePrimitiveColliders == true) UpdatePrims();
            if (timeSinceStartedComputing < minComputeTime)
            {
                timeSinceStartedComputing += Time.unscaledDeltaTime;
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

        #endregion Main

        #region API

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

        #endregion API

        #region ComputeMesh

        private ComputeMesh cm_job;
        private JobHandle cm_handle;
        private bool cm_isActive = false;
        private float timeSinceStartedComputing = 69420.0f;
        private float timeAtLastCompute = 0.0f;
        private float removeTimer = 0.0f;
        private int lastToRemoveCount = 0;

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

#if UNITY_EDITOR
            cm_job.mo.lToW = transform.localToWorldMatrix;
#endif

            float timeNow = Time.timeSinceLevelLoad;
            removeTimer += Mathf.Min(timeNow - timeAtLastCompute, minComputeTime * 1.2f);//Prevents removing all if got huge spike
            timeAtLastCompute = timeNow;
            float interval = voxRemoveInterval / Mathf.Max(1.0f, lastToRemoveCount * perVoxIntervalFactor);
            int remCount = Mathf.FloorToInt(removeTimer / interval);
            removeTimer -= interval * remCount;
            cm_job.mo.mayRemoveCount = remCount;

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
                                int vI = math.abs(voxI);
                                voxsToUpdate.Remove(vI);//If added and deleted before processed
                                if (voxIToCol.Remove(vI, out SphereCollider col) == false) continue;
                                col.gameObject.SetActive(false);
                                inactiveCols.Enqueue(col);
                                continue;
                            }

                            voxsToUpdate.Add(voxI);
                        }
                    }

                    if (Application.isPlaying == true)
                    {
                        ////Thanks to the way stuff is added to list, middle index is usually the rough center of all voxels that broke
                        //int centerVI = _mo.addedRemovedVoxsI[(_mo.addedRemovedVoxsI.Length - 1) / 2];
                        //if (centerVI >= 0) centerVI = _mo.addedRemovedVoxsI[^1];
                        //if (centerVI < 0)
                        //{
                        //    breakAudioProps.pos = transform.TransformPoint(_mo.VoxIndexToPosL(-centerVI));
                        //    breakAudio.Play(breakAudioProps);
                        //}

                        UpdateRigidbody();
                        //if (voxJobInitlized == true)
                        //    _voxsToSet_write.AddRange(_mo.addedRemovedVoxsI.AsArray());
                    }

                    _mo.addedRemovedVoxsI.Clear();
                }
            }
            else cm_job.mda.Dispose();

            if (_mo.voxsToCheck.Count > 0)
                chunksNeedsComputing = true;

            int remCount = _mo.voxsToRemove.Length;
            if (remCount > 0)
            {
                lastToRemoveCount = Mathf.Max(lastToRemoveCount, remCount);
                meshNeedsComputing = true;
            }
            else lastToRemoveCount = 0;

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
            public bool smoothNormals;
            public BreakMode breakMode;

            public void Execute()
            {
                if (chunksOnly == false)
                {
                    foreach (Action a in actions)
                    {
                        int voxI = mo.PosLToVoxIndex(wToL.MultiplyPoint3x4(a.posW + (0.5f * MeshGenGlobals.voxelSizeWorld * Vector3.one)));
                        Debug.DrawLine(a.posW, a.posW + Vector3.up, Color.red, 2.0f, false);
                        mo.SetVoxelsAround(voxI, a.radiusW, a.type == Action.Type.remove ? -1 : 1, a.radiusW < 2.01f);
                    }

                    mo.Meshify(mda[0], smoothNormals, out Bounds bl);
                    result.Value = new(bl);
                }

                if (breakMode == BreakMode.none)
                {
                    mo.voxsToCheck.Clear();
                    return;
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
            int vCountYZ = _mo.data.vCountYZ;
            int vCountZ = _mo.data.vCountZ;
            float3 voxSize = _mo.data.voxSize;

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
                    if ((vState & (1 << 1)) != 0) vI += vCountYZ;
                    if ((vState & (1 << 2)) != 0) vI += vCountZ;
                    if ((vState & (1 << 3)) != 0) vI += 1;

                    col.transform.localPosition = _mo.VoxIndexToPosL(vI);
                    continue;
                }

                col.transform.localPosition = _mo.VoxIndexToPosL(voxI)
                    + (voxSize * 0.5f);
            }

            if (loopCount >= maxPrimUpdatesPerFrame)
            {
                foreach (int voxI in voxsUpdated)
                {
                    voxsToUpdate.Remove(voxI);
                }
            }
            else voxsToUpdate.Clear();

            if (voxsToUpdate.Count == 0 && rb != null)//Freezing rb until all cols has been updated seems improve stability significantly
                rb.constraints = RigidbodyConstraints.None;
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
            transform.GetPositionAndRotation(out Vector3 pos, out Quaternion rot);
            newO.transform.SetPositionAndRotation(pos, rot);
            newO.layer = chunkBaseLayer;

            Rigidbody newRb = newO.GetComponent<Rigidbody>();
            newRb.interpolation = interpolation;
            newRb.constraints = RigidbodyConstraints.FreezeAll;
            //newRb.MoveRotation(rot);
            newRb.rotation = rot;//Setting both transform and rb seems to be the only reliable way to make sure it applies, wtf?
            //newRb.MovePosition(pos);
            newRb.position = pos;

            //if (this.rb != null)
            //{
            //    newRb.velocity = this.rb.velocity;
            //    newRb.angularVelocity = this.rb.angularVelocity;
            //    newRb.useGravity = this.rb.useGravity;
            //}

            MeshGenObject mgo = newO.GetComponent<MeshGenObject>();
            mgo.SetMO(new(mo));
            mgo.isChunk = true;
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
            mgo.smoothNormals = smoothNormals;

            tempDisableInit = false;
            mgo.Init();//~30% of creation time is this ~0.085ms (Editor)

            ignoreFrames = 2;
            mo.overlappingVoxCount = -1;
        }

        #endregion Chunks

        //#region Voxels
        //
        //private bool voxJobInitlized = false;
        //private ComputeVoxels cv_job;
        //private JobHandle cv_handle;
        //private bool voxsFlipped = false;
        //private NativeList<int> voxsToSetA;
        //private NativeList<int> voxsToSetB;
        //private NativeList<int> _voxsToSet_write => voxsFlipped == true ? voxsToSetB : voxsToSetA;
        //private NativeArray<int>.ReadOnly _voxsToSet_read => voxsFlipped == true
        //    ? voxsToSetA.AsReadOnly() : voxsToSetB.AsReadOnly();
        //
        //private void InitVoxels()
        //{
        //    if (voxJobInitlized == true) return;
        //    if (usePrimitiveColliders == true || Application.isPlaying == false) return;
        //    voxJobInitlized = true;
        //
        //    voxsToSetA = new(64, Allocator.Persistent);
        //    voxsToSetB = new(64, Allocator.Persistent);
        //    VoxGlobalHandler.OnScheduleWriteOperations += OnScheduleWriteOperations;
        //    VoxGlobalHandler.OnSetupVoxelSystem += OnSetupVoxelSystem;
        //    if (VoxGlobalHandler._hasBeenSetup == true) OnSetupVoxelSystem();
        //}
        //
        //private void DisposeVoxels()
        //{
        //    if (voxJobInitlized == false) return;
        //    voxJobInitlized = false;
        //
        //    VoxGlobalHandler.OnScheduleWriteOperations -= OnScheduleWriteOperations;
        //    VoxGlobalHandler.OnSetupVoxelSystem -= OnSetupVoxelSystem;
        //    cv_handle.Complete();
        //    voxsToSetA.Dispose();
        //    voxsToSetB.Dispose();
        //}
        //
        //private unsafe void OnSetupVoxelSystem()
        //{
        //    cv_job = new()
        //    {
        //        mod = _mo.data,
        //        lToW = transform.localToWorldMatrix,
        //        vWorld = VoxGlobalHandler._voxWorldNative_readonly,
        //        voxsCount = VoxGlobalHandler._voxGridCount_readWrite,
        //        voxsType = VoxGlobalHandler._voxGrid_readWrite,
        //        voxsTypeOld = VoxGlobalHandler._voxGridOld_readWrite,
        //    };
        //}
        //
        //private JobHandle OnScheduleWriteOperations(JobHandle dependOn)
        //{
        //    voxsFlipped = !voxsFlipped;
        //    cv_job.voxsToSet_read = _voxsToSet_read;
        //    _voxsToSet_write.Clear();
        //    cv_handle = cv_job.Schedule(dependOn);
        //    return cv_handle;
        //}
        //
        //[BurstCompile]
        //private unsafe struct ComputeVoxels : IJob
        //{
        //    public NativeArray<int>.ReadOnly voxsToSet_read;
        //    [NativeDisableUnsafePtrRestriction] public MarchingObject.Data mod;
        //    public Matrix4x4 lToW;
        //    public NativeReference<VoxWorld>.ReadOnly vWorld;
        //    [NativeDisableUnsafePtrRestriction] public byte* voxsCount;
        //    [NativeDisableUnsafePtrRestriction] public byte* voxsType;
        //    [NativeDisableUnsafePtrRestriction] public byte* voxsTypeOld;
        //
        //    public void Execute()
        //    {
        //        VoxWorld vWorld = this.vWorld.Value;
        //        int vwCountZ = vWorld.vCountZ;
        //        int vwCountZY = vWorld.vCountYZ;
        //        float worldMaxX = vWorld.vCountX * VoxGlobalSettings.voxelSizeWorld;
        //        float worldMaxY = vWorld.vCountY * VoxGlobalSettings.voxelSizeWorld;
        //        float worldMaxZ = vWorld.vCountZ * VoxGlobalSettings.voxelSizeWorld;
        //
        //        foreach (int mo_voxI in voxsToSet_read)
        //        {
        //            Vector3 voxPos = lToW.MultiplyPoint3x4(mod.VoxIndexToPosL(math.abs(mo_voxI)));
        //            if (voxPos.x < 0 || voxPos.y < 0 || voxPos.z < 0//Prevent out of bounds (Accurate, no warping)
        //                || voxPos.x > worldMaxX || voxPos.y > worldMaxY || voxPos.z > worldMaxZ) continue;
        //
        //            int wvIndex = (int)(voxPos.z * VoxGlobalSettings.voxelSizeWorldInv)
        //                + ((int)(voxPos.y * VoxGlobalSettings.voxelSizeWorldInv) * vwCountZ)
        //                + ((int)(voxPos.x * VoxGlobalSettings.voxelSizeWorldInv) * vwCountZY);
        //
        //            if (mo_voxI > 0)
        //            {
        //                VoxHelpBurst.AddVoxAtPos(wvIndex, VoxGlobalSettings.defualtType, voxsCount,
        //                    voxsType, voxsTypeOld, vwCountZ, vwCountZY);
        //            }
        //            else if (mo_voxI != 0)
        //            {
        //                VoxHelpBurst.RemoveVoxAtPos(wvIndex, VoxGlobalSettings.defualtType, voxsCount,
        //                    voxsType, voxsTypeOld, vwCountZ, vwCountZY);
        //            }
        //        }
        //    }
        //}
        //
        //#endregion Voxels
    }
}


