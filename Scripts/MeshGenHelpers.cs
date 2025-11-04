//Copyright 2025 David Westberg (MIT) https://github.com/Zombie1111/UnityMarchingCubes
using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using Unity.Collections;
using Unity.Collections.LowLevel.Unsafe;
using Unity.Mathematics;
using UnityEngine;
using UnityEngine.Rendering;

namespace zombGen
{
    public unsafe struct MarchingObject
    {
        public float* voxs;
        public volatile byte* voxsState;
        public float surface;
        public float iso;
        public bool hasPrimitives;

        public NativeList<uint> tris;
        public NativeList<float3> vers;
        public NativeList<float3> nors;
        public NativeReference<int> solidVoxCount;
        public Data data;

        public struct Data
        {
            /// <summary>
            /// Lenght of voxs, -1 if invalid
            /// </summary>
            public int vCountXYZ;
            public int vCountYZ;
            public int vCountZ;
            public int vCountY;
            public int vCountX;
            public float3 start;
            public float3 voxSize;

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public readonly float3 VoxIndexToPosL(int voxI)
            {
                int remainderAfterZ = voxI % vCountYZ;
                return new float3((voxI / vCountYZ) * voxSize.x,
                    (remainderAfterZ / vCountZ) * voxSize.y,
                    (remainderAfterZ % vCountZ) * voxSize.z) + start;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public readonly int PosLToVoxIndex(float3 posL, int voxMargin = 1)
            {
                posL -= start;

                float3 snappedPos = new(
                    math.clamp(posL.x, voxSize.x * 1.5f * voxMargin, (vCountX * voxSize.x) - (voxSize.x * 1.5f * voxMargin)),
                    math.clamp(posL.y, voxSize.y * 1.5f * voxMargin, (vCountY * voxSize.y) - (voxSize.y * 1.5f * voxMargin)),
                    math.clamp(posL.z, voxSize.z * 1.5f * voxMargin, (vCountZ * voxSize.z) - (voxSize.z * 1.5f * voxMargin))
                );

                int vI = (int)(snappedPos.z / voxSize.z)
                    + ((int)(snappedPos.y / voxSize.y) * vCountZ)
                    + ((int)(snappedPos.x / voxSize.x) * vCountYZ);
                return vI;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public readonly bool IsVoxIndexInBounds(int voxI)
            {
                int remainderAfterZ = voxI % vCountYZ;
                int x = voxI / vCountYZ;
                int y = remainderAfterZ / vCountZ;
                int z = remainderAfterZ % vCountZ;
                return x >= 0 && x < vCountX
                    && y >= 0 && y < vCountY
                    && z >= 0 && z < vCountZ;
            }

            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            public readonly bool IsVoxIndexInBounds(int voxI, int margin)
            {
                int remainderAfterZ = voxI % vCountYZ;
                int x = voxI / vCountYZ;
                int y = remainderAfterZ / vCountZ;
                int z = remainderAfterZ % vCountZ;
                return x >= margin && x < vCountX - margin
                    && y >= margin && y < vCountY - margin
                    && z >= margin && z < vCountZ - margin;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly float3 VoxIndexToPosL(int voxI)
        {
            return data.VoxIndexToPosL(voxI);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly int PosLToVoxIndex(float3 posL, int voxMargin = 1)
        {
            return data.PosLToVoxIndex(posL, voxMargin);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool IsVoxIndexInBounds(int voxI)
        {
            return data.IsVoxIndexInBounds(voxI);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool IsVoxIndexInBounds(int voxI, int margin)
        {
            return data.IsVoxIndexInBounds(voxI, margin);
        }

        /// <summary>
        /// Voxel indexs of removed or created voxels (Negative if removed)
        /// DONT forget to clear
        /// </summary>
        public NativeList<int> addedRemovedVoxsI;//voxel 0 should never be able to be solid so no need for +1

        #region Voxelization

        /// <summary>
        /// DO NOT ACCESS this after disposing (Returns new cleared VoxObject)
        /// </summary>
        public void Dispose()
        {
            if (data.vCountXYZ < 0) return;
            data.vCountXYZ = -1;

            UnsafeUtility.Free(voxs, Allocator.Persistent);
            UnsafeUtility.Free(voxsState, Allocator.Persistent);
            tris.Dispose();
            vers.Dispose();
            nors.Dispose();
            addedRemovedVoxsI.Dispose();
            newMO.Value.Dispose();
            newMO.Dispose();
            voxsToCheck.Dispose();
            solidVoxCount.Dispose();
            voxsToRemove.Dispose();
        }

        /// <summary>
        /// Voxelizes the given collider. (Mainthread only)
        /// DONT FORGET TO CALL this.Dispose() in OnDestroy or similar
        /// </summary>
        public MarchingObject(Collider col, Collider kin = null, bool tryFillConvaveInteriors = true, bool tryStripOuterLayer = false,
            float surface = 0.0f, float isoLevel = 1.0f, bool hasPrimitives = false,
            MeshGenObject.BreakMode breakMode = MeshGenObject.BreakMode.rigidbody)
        {
#if UNITY_EDITOR
            lToW = col.transform.localToWorldMatrix;
#endif

            Quaternion ogRot = col.transform.rotation;
            col.transform.rotation = Quaternion.identity;
            col.enabled = false;
            col.enabled = true;

            Vector3 voxelSize = MeshGenGlobals.voxelSizeWorld * Vector3.one;
            Bounds colBounds = col.bounds;
            Matrix4x4 colWToL = col.transform.worldToLocalMatrix;
            Vector3 inverseScale = col.transform.lossyScale;
            inverseScale = new Vector3(1.0f / inverseScale.x, 1.0f / inverseScale.y, 1.0f / inverseScale.z);

            Vector3 bSize = colBounds.size + (voxelSize * 2.0f);
            int vCountX = (int)Math.Ceiling(bSize.x / voxelSize.x);
            int vCountY = (int)Math.Ceiling(bSize.y / voxelSize.y);
            int vCountZ = (int)Math.Ceiling(bSize.z / voxelSize.z);
            int vCountYZ = vCountY * vCountZ;
            bool hasAnyKin = false;

            Vector3 bStart = colBounds.min - voxelSize;
            Vector3 bMax = colBounds.max;
            int ogLayer = col.gameObject.layer;
            int kinOgLayer = kin == null ? 0 : kin.gameObject.layer;
            bool kinActive = kin != null && kin.gameObject.activeSelf;
            bool kinEnabled = kin != null && kin.enabled;
            bool kinTrigger = kin != null && kin.isTrigger;
            col.gameObject.layer = MeshGenGlobals.voxelizeTempLayer;
            if (kin != null)
            {
                kin.gameObject.layer = MeshGenGlobals.voxelizeTempLayer;
                kin.isTrigger = false;
                kin.gameObject.SetActive(true);
                kin.enabled = true;
            }
            if (ogLayer == MeshGenGlobals.voxelizeTempLayer) Debug.LogWarning(col.transform.name + " already had voxelizeTempLayer layer, not recommended");
            LayerMask layerMask = 1 << MeshGenGlobals.voxelizeTempLayer;
            PhysicsScene colPhyScene = col.gameObject.scene.GetPhysicsScene();

            //Overlap tests
            Vector3 voxelSizeHalf = voxelSize * 0.501f;
            Vector3 voxelSizeHalfReal = voxelSize * 0.5f;
            int vCountXYZ = vCountX * vCountY * vCountZ;
            float* voxs = null;
            byte* voxsState = null;
            int overlappingVoxCount = 0;
            int solidVoxCount = 0;

            if (vCountXYZ < 9)
            {
                vCountXYZ = -1;
                Debug.LogWarning(col.transform.name + " wont contain any voxels! Is the object too small?");
                goto SkipCreation;
            }

            if (vCountXYZ > MeshGenGlobals.maxVoxelsInExtent
                * MeshGenGlobals.maxVoxelsInExtent
                * MeshGenGlobals.maxVoxelsInExtent)
            {
                vCountXYZ = -1;
                Debug.Log(col.transform.name + " cant be voxelized because its bounds are too large, bounds volume: " + vCountXYZ);
                goto SkipCreation;
            }

            var commands = new NativeArray<OverlapBoxCommand>(vCountXYZ, Allocator.TempJob);
            var results = new NativeArray<ColliderHit>(vCountXYZ * MeshGenGlobals.voxelizeMaxHits, Allocator.TempJob);
            QueryParameters qParams = new() { layerMask = layerMask, hitBackfaces = true, hitMultipleFaces = true, hitTriggers = QueryTriggerInteraction.Ignore };

            Parallel.For(0, vCountX, x =>
            {
                int index = x * vCountYZ;

                for (int y = 0; y < vCountY; y++)
                {
                    for (int z = 0; z < vCountZ; z++)
                    {
                        commands[index] = new OverlapBoxCommand((new Vector3(x, y, z) * MeshGenGlobals.voxelSizeWorld) + bStart, voxelSizeHalf,
                            Quaternion.identity, qParams);
                        index++;
                    }
                }
            });

            OverlapBoxCommand.ScheduleBatch(commands, results, 1, MeshGenGlobals.voxelizeMaxHits).Complete();
            commands.Dispose();

            //Check result
            int colId = col.GetInstanceID();
            int kinId = kin == null ? 0 : kin.GetInstanceID();

            voxs = (float*)UnsafeUtility.Malloc(vCountXYZ * sizeof(float),
                UnsafeUtility.AlignOf<float>(), Allocator.Persistent);
            voxsState = (byte*)UnsafeUtility.Malloc(vCountXYZ * sizeof(byte),
                UnsafeUtility.AlignOf<byte>(), Allocator.Persistent);
            UnsafeUtility.MemClear(voxsState, vCountXYZ * sizeof(byte));

            Parallel.For(0, vCountXYZ, vI =>
            {
                int maxRI = (vI * MeshGenGlobals.voxelizeMaxHits) + MeshGenGlobals.voxelizeMaxHits;
                for (int rI = vI * MeshGenGlobals.voxelizeMaxHits; rI < maxRI; rI++)
                {
                    int id = results[rI].instanceID;
                    if (id == kinId)
                    {
                        hasAnyKin = true;
                        voxsState[vI] = MeshGenGlobals.kinVoxFlag;
                    }
                    if (id != colId) continue;


                    voxs[vI] = isoLevel;
                    overlappingVoxCount++;
                    return;
                }

                voxs[vI] = -1.0f;
            });

            solidVoxCount = overlappingVoxCount;
            results.Dispose();

            //Floodfill
            if (tryFillConvaveInteriors == true && col is MeshCollider mc && mc.convex == false && mc.sharedMesh != null)
            {
                Mesh m = mc.sharedMesh;
                int[] tris = m.triangles;
                Vector3[] vers = m.vertices;
                int triCount = tris.Length;
                Matrix4x4 lToW = mc.transform.localToWorldMatrix;

                for (int i = 0; i < triCount; i += 3)//We can skip every other, odds of this missing a hole is basically 0
                {//No, skipping every other will probably decrease performance since it will result in more iterations
                    Vector3 v0 = vers[tris[i]];
                    Vector3 v1 = vers[tris[i + 1]];
                    Vector3 v2 = vers[tris[i + 2]];
                    Vector3 nor = Vector3.Normalize(lToW.MultiplyVector(Vector3.Cross(v1 - v0, v2 - v0)));
                    Vector3 center = lToW.MultiplyPoint3x4((v0 + v1 + v2) / 3f) + (nor * MeshGenGlobals.voxelSizeWorld);

                    if (center.x < bStart.x || center.x > bMax.x
                        || center.y < bStart.y || center.y > bMax.y
                        || center.z < bStart.z || center.z > bMax.z)
                        continue;

                    center += voxelSizeHalfReal;
                    center -= bStart;

                    int vI = (int)(center.z / MeshGenGlobals.voxelSizeWorld)
                        + ((int)(center.y / MeshGenGlobals.voxelSizeWorld) * vCountZ)
                        + ((int)(center.x / MeshGenGlobals.voxelSizeWorld) * vCountYZ);
                    if (voxs[vI] != -1.0f) continue;
                    voxs[vI] = 2.0f;
                }

                bool spreadedAny = true;
                int vMaxX = vCountX - 1;
                int vMaxY = vCountY - 1;
                int vMaxZ = vCountZ - 1;
                int[] offsets = new int[6]
                {
                    1, -1,
                    vCountZ, -vCountZ,
                    vCountYZ, -vCountYZ
                };

                while (spreadedAny == true)
                {
                    spreadedAny = false;

                    for (int vI = 0; vI < vCountXYZ; vI++)
                    {
                        if (voxs[vI] != 2.0f) continue;
                        voxs[vI] = 3.0f;

                        for (int i = 0; i < 6; i++)
                        {
                            int nVI = vI + offsets[i];
                            int remainderAfterZ = nVI % vCountYZ;
                            int x = nVI / vCountYZ;
                            int y = remainderAfterZ / vCountZ;
                            int z = remainderAfterZ % vCountZ;
                            if (x < 0 || x > vMaxX//Bounds check
                                || y < 0 || y > vMaxY
                                || z < 0 || z > vMaxZ) continue;

                            if (voxs[nVI] != -1.0f) continue;
                            voxs[nVI] = 2.0f;
                            spreadedAny = true;
                        }
                    }
                }

                if (tryStripOuterLayer == true)
                {
                    for (int vI = 0; vI < vCountXYZ; vI++)
                    {
                        if (voxs[vI] != isoLevel) continue;
                        voxs[vI] = 3.0f;
                    }
                }

                for (int vI = 0; vI < vCountXYZ; vI++)
                {
                    if (voxs[vI] == -1.0f)
                    {
                        solidVoxCount++;
                        voxs[vI] = isoLevel;
                    }
                    else if (voxs[vI] == 3.0f) voxs[vI] = -1.0f;
                }
            }

        SkipCreation:;

            this.voxs = voxs;
            this.voxsState = voxsState;
            data = new()
            {
                start = colWToL.MultiplyPoint3x4(bStart),
                voxSize = inverseScale * MeshGenGlobals.voxelSizeWorld,
                vCountXYZ = vCountXYZ,
                vCountYZ = vCountYZ,
                vCountZ = vCountZ,
                vCountY = vCountY,
                vCountX = vCountX,
            };

            this.surface = surface;
            this.hasAnyKin = hasAnyKin;
            iso = isoLevel;
            this.hasPrimitives = hasPrimitives;
            mayRemoveCount = 0;
            removeFloatingVoxs = breakMode == MeshGenObject.BreakMode.remove;
            voxsToRemove = new(removeFloatingVoxs == false ? 1 : 128, Allocator.Persistent);

            addedRemovedVoxsI = new(overlappingVoxCount, Allocator.Persistent);
            newMO = new(new(removeFloatingVoxs == true ? 9 : vCountXYZ), Allocator.Persistent);
            voxsToCheck = new(32, Allocator.Persistent);
            this.solidVoxCount = new(solidVoxCount, Allocator.Persistent);

            int initialVerCount = overlappingVoxCount * 5 * 3;
            tris = new(initialVerCount, Allocator.Persistent);
            vers = new(initialVerCount, Allocator.Persistent);
            nors = new(initialVerCount, Allocator.Persistent);

            col.gameObject.layer = ogLayer;
            if (kin != null)
            {
                kin.gameObject.layer = kinOgLayer;
                kin.gameObject.SetActive(kinActive);
                kin.enabled = kinEnabled;
                kin.isTrigger = kinTrigger;
            }
            col.transform.rotation = ogRot;
        }
        #endregion Voxelization

        #region Meshification

        public int mayRemoveCount;

        public void Meshify(Mesh.MeshData md, bool smoothNormals, out Bounds boundsLocal)
        {
            int vCountZ = data.vCountZ;
            int vCountX = data.vCountX;
            int vCountY = data.vCountY;
            int vCountYZ = data.vCountYZ;
            int vCountXYZ = data.vCountXYZ;
            bool hadStuffToCheck = voxsToCheck.Count > 0;

            //Remove voxs
            int remCount = math.min(voxsToRemove.Length, mayRemoveCount);
            if (remCount > 0)
            {
                for (int i = 0; i < remCount; i++)
                {
                    int vI = voxsToRemove[i];
                    voxs[vI] = -1.0f;
                    voxsState[vI] = (byte)(voxsState[vI] & ~MeshGenGlobals.delVoxFlag);
                }

                voxsToRemove.RemoveRange(0, remCount);
            }

            //Alloc
            float* cubeValues = (float*)UnsafeUtility.Malloc(8 * sizeof(float),
                UnsafeUtility.AlignOf<float>(), Allocator.Temp);
            float3* edgeVertex = (float3*)UnsafeUtility.Malloc(12 * 3 * sizeof(float),
                UnsafeUtility.AlignOf<float3>(), Allocator.Temp);

            //We could persistent allocate these in a readonly struct and pass along that, worth it?
            int* offsets = MeshGenLookup.GetOffsets(vCountZ, vCountYZ);
            int* offsetsAll = MeshGenLookup.GetOffsetsAll(vCountZ, vCountYZ);
            int* cubeEdgeFlags = MeshGenLookup.GetCubeEdgeFlags();
            int2* edgeConnections = MeshGenLookup.GetEdgeConnections();
            float3* edgeDirections = MeshGenLookup.GetEdgeDirections();
            float3* vertexOffsets = MeshGenLookup.GetVertexOffsets();
            int* triangleConnectionTables = MeshGenLookup.GetTriangleConnectionTables();
            UnsafeHashMap<int3, uint> verNorMap = new(Mathf.Max(solidVoxCount.Value * 3, 16), Allocator.Temp);

            var _vers = vers;
            var _tris = tris;
            var _nors = nors;

            float3 voxSize = data.voxSize;
            float3 minPos = new(69420.0f, 69420.0f, 69420.0f);
            float3 maxPos = -minPos;
            bool hadMesh = vers.Length > 2;
            vers.Clear();
            tris.Clear();
            nors.Clear();
            uint verLenght = 0;

            for (int vI = 0; vI < vCountXYZ; vI++)
            {
                float3 posL = VoxIndexToPosL(vI);
                int cubeIndex = 0;
                int vMaxX = vCountX - 1;
                int vMaxY = vCountY - 1;
                int vMaxZ = vCountZ - 1;

                for (int i = 0; i < 8; i++)
                {
                    int nVI = vI + offsets[i];
                    if (IsVoxIndexInBounds(nVI) == false)
                    {//Bounds check
                        cubeValues[i] = -1.0f;
                        continue;
                    }

                    cubeValues[i] = voxs[nVI];
                    if (IsVoxIndexSolid(nVI) == true)
                        cubeIndex |= 1 << i;
                }

                int edgeFlag = cubeEdgeFlags[cubeIndex];
                if (edgeFlag == 0)
                {
                    if ((voxsState[vI] & ~MeshGenGlobals.kinDelVoxFlag) == 0) continue;
                    voxsState[vI] &= MeshGenGlobals.kinDelVoxFlag;
                    addedRemovedVoxsI.Add(-vI);

                    for (int i = 0; i < 14; i++)
                    {
                        int vI3 = vI + offsetsAll[i];
                        if (IsVoxIndexInBounds(vI3) == false) continue;
                        if ((voxsState[vI3] & MeshGenGlobals.delVoxFlag) == 0)
                            voxsToCheck.Add(vI3);
                    }

                    continue;
                }

                for (int i = 0; i < 12; i++)
                {
                    //if there is an intersection on this edge
                    if ((edgeFlag & (1 << i)) == 0) continue;
                    int2 edgeCon = edgeConnections[i];
                    float offset = GetOffset(cubeValues[edgeCon.x],
                        cubeValues[edgeCon.y]);

                    edgeVertex[i] = posL + ((vertexOffsets[edgeCon.x] + offset * edgeDirections[i]) * voxSize);
                }

                float3 totNor = float3.zero;
                int cubeIndex_tri = cubeIndex * 16;
                for (int i = 0; i < 16; i += 3)
                {
                    if (triangleConnectionTables[cubeIndex_tri + i] < 0) break;

                    float3 ev0 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i]];
                    float3 ev1 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i + 1]];
                    float3 ev2 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i + 2]];
                    //float3 nor = math.normalize(math.cross(ev1 - ev2, ev0 - ev2));
                    float3 nor = math.cross(ev1 - ev2, ev0 - ev2);//normalize seems to have little to no affect
                    int3 norKey = smoothNormals == false ? (int3)(nor * 10000) : int3.zero;
                    totNor += nor;
                    Add(ev2);
                    Add(ev1);
                    Add(ev0);

                    [MethodImpl(MethodImplOptions.AggressiveInlining)]
                    void Add(float3 ver)
                    {
                        int3 key = norKey + (int3)(ver * 10000);
                        if (verNorMap.TryGetValue(key, out uint idx) == false)
                        {
                            idx = verLenght++;
                            verNorMap.Add(key, idx);

                            _vers.Add(ver);
                            _nors.Add(nor);
                            minPos = math.min(minPos, math.min(ev0, math.min(ev1, ev2)));
                            maxPos = math.max(maxPos, math.max(ev0, math.max(ev1, ev2)));
                        }
                        else
                        {
                            _nors[(int)idx] += nor;
                        }

                        _tris.Add(idx);
                    }
                }

                //Could be precomputed in lookup table based on cubeIndex
                byte prevState = voxsState[vI];
                byte newState = (byte)((1 << 0) | (prevState & MeshGenGlobals.kinDelVoxFlag));
                if (totNor.x < 0.0f) newState |= 1 << 1;
                if (totNor.y < 0.0f) newState |= 1 << 2;
                if (totNor.z < 0.0f) newState |= 1 << 3;

                if (prevState == newState) continue;
                voxsState[vI] = newState;
                if ((prevState & ~MeshGenGlobals.kinDelVoxFlag) == 0
                    || hasPrimitives == true) addedRemovedVoxsI.Add(vI);

                //if ((prevState & ~MeshGenGlobals.kinDelVoxFlag) == 0 && hadMesh == true)
                if (hadMesh == true)
                {
                    int vI2 = vI;
                    if ((newState & (1 << 1)) != 0) vI2 += vCountYZ;
                    if ((newState & (1 << 2)) != 0) vI2 += vCountZ;
                    if ((newState & (1 << 3)) != 0) vI2 += 1;

                    if ((voxsState[vI2] & MeshGenGlobals.delVoxFlag) == 0)
                        voxsToCheck.Add(vI2);

                    //for (int i = 0; i < 14; i++)
                    //{
                    //    int vI3 = vI2 + offsetsAll[i];
                    //    if (IsVoxIndexInBounds(vI3) == false) continue;
                    //    if ((voxsState[vI3] & MeshGenGlobals.delVoxFlag) == 0)
                    //        voxsToCheck.Add(vI3);
                    //}
                }
            }

            //for (int i = 0; i < verLenght; i++)
            //{//Seems to work without normalizing
            //    _nors[i] = math.normalize(_nors[i]);
            //}

            //Write result to mesh
            boundsLocal = new((minPos + maxPos) / 2, maxPos - minPos);

            NativeArray<VertexAttributeDescriptor> layout = new(2, Allocator.Temp);
            layout[0] = new VertexAttributeDescriptor(VertexAttribute.Position, VertexAttributeFormat.Float32, 3, 0);
            layout[1] = new VertexAttributeDescriptor(VertexAttribute.Normal, VertexAttributeFormat.Float32, 3, 1);
            //layout[2] = new VertexAttributeDescriptor(VertexAttribute.Tangent, VertexAttributeFormat.Float32, 4, 2);
            //layout[3] = new VertexAttributeDescriptor(VertexAttribute.TexCoord0, VertexAttributeFormat.Float32, 2, 3);
            md.SetVertexBufferParams(vers.Length, layout);
            md.SetIndexBufferParams(tris.Length, IndexFormat.UInt32);
            md.subMeshCount = 1;
            UnsafeUtility.MemCpy(md.GetVertexData<Vector3>(0).GetUnsafePtr(),
                vers.GetUnsafeReadOnlyPtr(), vers.Length * sizeof(float3));
            UnsafeUtility.MemCpy(md.GetVertexData<Vector3>(1).GetUnsafePtr(),
                nors.GetUnsafeReadOnlyPtr(), nors.Length * sizeof(float3));
            UnsafeUtility.MemCpy(md.GetIndexData<uint>().GetUnsafePtr(),
                tris.GetUnsafeReadOnlyPtr(), tris.Length * sizeof(uint));
            md.SetSubMesh(0, new(0, tris.Length, MeshTopology.Triangles));

            //Dispose
            layout.Dispose();
            UnsafeUtility.Free(cubeValues, Allocator.Temp);
            UnsafeUtility.Free(edgeVertex, Allocator.Temp);
            UnsafeUtility.Free(offsets, Allocator.Temp);
            UnsafeUtility.Free(offsetsAll, Allocator.Temp);
            UnsafeUtility.Free(cubeEdgeFlags, Allocator.Temp);
            UnsafeUtility.Free(edgeConnections, Allocator.Temp);
            UnsafeUtility.Free(edgeDirections, Allocator.Temp);
            UnsafeUtility.Free(vertexOffsets, Allocator.Temp);
            UnsafeUtility.Free(triangleConnectionTables, Allocator.Temp);
            verNorMap.Dispose();

            ////Temp fix for floating shit
            ////if (hadMesh == true && voxsToCheck.Count == 0)
            //if (hadMesh == true && hadStuffToCheck == false)
            //{
            //    for (int i = 0; i < vCountXYZ; i++)
            //    {
            //        if (IsVoxIndexSolid(i) == false
            //            || (voxsState[i] & MeshGenGlobals.delVoxFlag) != 0) continue;
            //        voxsToCheck.Add(i);
            //    }
            //}
        }

        #endregion Meshification

        #region Base Helpers



        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private readonly float GetOffset(float v1, float v2)
        {
            float delta = v2 - v1;
            return (delta == 0.0f) ? surface : (surface - v1) / delta;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool IsVoxIndexSolid(int voxI)
        {
            return voxs[voxI] >= surface;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly bool IsVoxIndexKinematic(int voxI)
        {
            return (voxsState[voxI] & MeshGenGlobals.kinVoxFlag) != 0;
        }

        #endregion Base Helpers

        #region Methods

        public void SetVoxelsAround(int voxI, float radius, float newValue = -1.0f, bool setClosest = false)
        {
            radius /= 2.0f;//Looks more correct, idk
            float3 voxSize = data.voxSize;
            int vCountZ = data.vCountZ;
            int vCountYZ = data.vCountYZ;
            int maxX = Mathf.CeilToInt(radius / MeshGenGlobals.voxelSizeWorld);
            int maxY = Mathf.CeilToInt(radius / MeshGenGlobals.voxelSizeWorld);
            int maxZ = Mathf.CeilToInt(radius / MeshGenGlobals.voxelSizeWorld);
            float radiusSQ = radius * radius;
            bool newSolid = newValue >= surface;
            int margin = newSolid == true ? 1 : 0;
            int changeCount = 0;
            //Margin when creating as marching cubes requires 1 row of empty voxels to work properly

            float closestDis = float.MaxValue;

            for (int z = -maxZ; z <= maxZ; z++)
            {
                for (int y = -maxY; y <= maxY; y++)
                {
                    for (int x = -maxX; x <= maxX; x++)
                    {
                        float dxw = x * MeshGenGlobals.voxelSizeWorld;
                        float dyw = y * MeshGenGlobals.voxelSizeWorld;
                        float dzw = z * MeshGenGlobals.voxelSizeWorld;
                        float disSQ = dxw * dxw + dyw * dyw + dzw * dzw;
                        if (disSQ > radiusSQ) continue;

                        int vI = voxI + z + (vCountZ * y) + (vCountYZ * x);
                        if (IsVoxIndexInBounds(vI, margin) == false) continue;
                        if (setClosest == true && closestDis < disSQ) continue;
                        if (IsVoxIndexSolid(vI) != newSolid)
                        {
                            changeCount++;
                            closestDis = disSQ * 1.1f;
                        }

                        voxs[vI] = newValue;
                    }
                }
            }

            if (newSolid == false && changeCount > 0)
                solidVoxCount.Value -= changeCount;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void RemoveVoxelsAround(float3 posL, float radius)
        {
            SetVoxelsAround(PosLToVoxIndex(posL), radius, -1.0f);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public void CreateVoxelsAround(float3 posL, float radius)
        {
            SetVoxelsAround(PosLToVoxIndex(posL), radius, 1.0f);
        }
        #endregion Methods

        #region Physics

        public NativeHashSet<int> voxsToCheck;
        public NativeReference<Unmanaged> newMO;
        public bool hasAnyKin;
        public bool removeFloatingVoxs;
        public NativeList<int> voxsToRemove;
#if UNITY_EDITOR
        /// <summary>
        /// Editor only, used for debug
        /// </summary>
        public Matrix4x4 lToW;
#endif

        public void TryGetChunk()
        {
        TryAgain:;
            if (voxsToCheck.Count == 0) return;
            int voxI = 0;
            foreach (int vI in voxsToCheck)
            {
                voxI = vI;//Really the best way?
                break;
            }

            //int voxI = voxsToCheck[^1];
            //voxsToCheck.RemoveAt(voxsToCheck.Length - 1);
            voxsToCheck.Remove(voxI);

            if (IsVoxIndexSolid(voxI) == false)
            {
                goto TryAgain;
            }

            int maxSpread = MeshGenGlobals.maxChunkCheckRadius - 2;
            int maxLoops = maxSpread * maxSpread * maxSpread;
            NativeHashSet<int> voxsSearched = new(maxLoops, Allocator.Temp);
            NativeQueue<int> toSearch = new(Allocator.Temp);
            int* offsets = MeshGenLookup.GetOffsetsStraight(data.vCountZ, data.vCountYZ);
            int ogToRemoveLen = voxsToRemove.Length;
            int loopCount = 0;
            bool isKin = false;

            toSearch.Enqueue(voxI);
            voxsSearched.Add(voxI);
            voxsToRemove.Add(voxI);

            while (loopCount++ < maxLoops && toSearch.TryDequeue(out int vI) == true)
            {
                for (int i = 0; i < 6; i++)
                {
                    int nVI = vI + offsets[i];
                    if (IsVoxIndexInBounds(nVI) == false) continue;
                    if (IsVoxIndexSolid(nVI) == false) continue;
                    if (voxsSearched.Add(nVI) == false) continue;

                    voxsToCheck.Remove(nVI);
                    voxsToRemove.Add(nVI);
                    //if (loopCount < 256)//As stuff likely changed in the same area, only check in the first few iterations 
                    //{
                    //    int toCheckI = voxsToCheck.IndexOf(nVI);
                    //    if (toCheckI >= 0) voxsToCheck.RemoveAtSwapBack(toCheckI);
                    //}

                    if (IsVoxIndexKinematic(nVI) == true)
                    {
                        isKin = true;
                        continue;
                    }

                    Debug.DrawLine(lToW.MultiplyPoint3x4(VoxIndexToPosL(vI)), lToW.MultiplyPoint3x4(VoxIndexToPosL(nVI)), Color.red, 0.0f, true);
                    toSearch.Enqueue(nVI);
                }
            }

            int searchedCount = voxsSearched.Count;
            toSearch.Dispose();
            UnsafeUtility.Free(offsets, Allocator.Temp);
            if (loopCount >= maxLoops || isKin == true || (hasAnyKin == false
                && searchedCount > solidVoxCount.Value - searchedCount))
            {
                voxsToRemove.Length = ogToRemoveLen;
                voxsSearched.Dispose();
                return;
            }

            //We have a disconnected chunk
            if (removeFloatingVoxs == true)
            {
                for (int i = voxsToRemove.Length - 1; i >= ogToRemoveLen; i--)
                {
                    voxsState[voxsToRemove[i]] |= MeshGenGlobals.delVoxFlag;
                }

                voxsSearched.Dispose();
                return;
            }

            Unmanaged newMO = this.newMO.Value;
            newMO.SetFrom(ref this, voxsSearched);
            this.newMO.Value = newMO;
            voxsSearched.Dispose();
        }

        public unsafe struct Unmanaged
        {
            public float* voxs;
            public byte* voxsState;
            public float3 start;
            public float3 voxSize;
            public float surface;
            public float iso;

            /// <summary>
            /// Lenght of voxs, -1 if invalid
            /// </summary>
            public int maxVCountXYZ;
            public int vCountXYZ;
            public int vCountYZ;
            public int vCountZ;
            public int vCountY;
            public int vCountX;

            public int overlappingVoxCount;
            public int solidVoxCount;

            public void Dispose()
            {
                if (maxVCountXYZ < 0) return;
                maxVCountXYZ = -1;
                UnsafeUtility.Free(voxs, Allocator.Persistent);
                UnsafeUtility.Free(voxsState, Allocator.Persistent);
            }

            public Unmanaged(int vCountXYZ)
            {
                maxVCountXYZ = Mathf.Min(vCountXYZ, MeshGenGlobals.maxChunkCheckRadius
                    * MeshGenGlobals.maxChunkCheckRadius
                    * MeshGenGlobals.maxChunkCheckRadius);
                if (maxVCountXYZ < 9) throw new Exception("Bad input " + vCountXYZ);

                this.vCountXYZ = -1;
                vCountYZ = -1;
                vCountZ = -1;
                vCountY = -1;
                vCountX = -1;
                overlappingVoxCount = -1;
                solidVoxCount = -1;

                start = float3.zero;
                voxSize = float3.zero;
                surface = 0.0f;
                iso = 1.0f;

                voxs = (float*)UnsafeUtility.Malloc(maxVCountXYZ * sizeof(float),
                    UnsafeUtility.AlignOf<float>(), Allocator.Persistent);
                voxsState = (byte*)UnsafeUtility.Malloc(maxVCountXYZ * sizeof(byte),
                    UnsafeUtility.AlignOf<byte>(), Allocator.Persistent);
            }

            internal void SetFrom(ref MarchingObject mo, NativeHashSet<int> voxsSearched)
            {
                //Get bounds
                NativeArray<float3> voxsPos = new(voxsSearched.Count, Allocator.Temp, NativeArrayOptions.UninitializedMemory);
                float3 bMin = new(69420.0f, 69420.0f, 69420.0f);
                float3 bMax = -bMin;
                int i = 0;

                foreach (int vI in voxsSearched)
                {
                    float3 posL = mo.VoxIndexToPosL(vI);
                    bMin = math.min(posL, bMin);
                    bMax = math.max(posL, bMax);
                    voxsPos[i++] = posL;
                }

                overlappingVoxCount = voxsPos.Length;
                solidVoxCount = voxsPos.Length;
                float3 voxelSize = mo.data.voxSize;
                float3 bSize = (bMax - bMin) + (voxelSize * (solidVoxCount > 1 ? 2 : 3));
                float3 bStart = bMin - voxelSize;

                vCountX = (int)Math.Ceiling(bSize.x / voxelSize.x);
                vCountY = (int)Math.Ceiling(bSize.y / voxelSize.y);
                vCountZ = (int)Math.Ceiling(bSize.z / voxelSize.z);
                vCountYZ = vCountY * vCountZ;
                vCountXYZ = vCountX * vCountY * vCountZ;
                voxSize = mo.data.voxSize;
                surface = mo.surface;
                iso = mo.iso;
                start = bStart;

                if (vCountXYZ < 9 || vCountXYZ > maxVCountXYZ)
                {
                    overlappingVoxCount = -1;
                    solidVoxCount = -1;
                    return;
                }

                mo.solidVoxCount.Value -= voxsPos.Length;
                foreach (int vI in voxsSearched)
                {
                    mo.voxs[vI] = -1.0f;//Remove voxel from source
                }

                float fill = -1.0f;
                UnsafeUtility.MemCpyReplicate(voxs, &fill, sizeof(float), vCountXYZ);
                UnsafeUtility.MemClear(voxsState, vCountXYZ * sizeof(byte));

                float3 offset = bStart - mo.data.start;
                i = 0;

                foreach (float3 vPos in voxsPos)
                {
                    //if (i++ > voxsPos.Length - 20) break;

                    //Debug.DrawRay(mo.lToW.MultiplyPoint3x4(vPos), Vector3.up * voxSize, Color.magenta, 5.0f, true);
                    float3 snappedPos = (vPos - mo.data.start) - offset;
                    int vI = (int)Math.Round(snappedPos.z / voxSize.z)
                        + ((int)Math.Round(snappedPos.y / voxSize.y) * vCountZ)
                        + ((int)Math.Round(snappedPos.x / voxSize.x) * vCountYZ);

                    voxs[vI] = mo.iso;
                }

                voxsPos.Dispose();
            }
        }

        public MarchingObject(in Unmanaged mo)
        {
            if (mo.overlappingVoxCount <= 0 || mo.solidVoxCount <= 0)
                throw new Exception("Bad Unmanaged overlappingVoxCount or solidVoxCount");

            data = new()
            {
                start = mo.start,
                voxSize = mo.voxSize,
                vCountXYZ = mo.vCountXYZ,
                vCountYZ = mo.vCountYZ,
                vCountZ = mo.vCountZ,
                vCountY = mo.vCountY,
                vCountX = mo.vCountX,
            };

            surface = mo.surface;
            iso = mo.iso;
            hasAnyKin = false;
            hasPrimitives = true;
            removeFloatingVoxs = false;
            voxsToRemove = new(1, Allocator.Persistent);
            mayRemoveCount = 0;

            int vCountXYZ = mo.vCountXYZ;
            voxs = (float*)UnsafeUtility.Malloc(sizeof(float) * vCountXYZ, UnsafeUtility.AlignOf<float>(), Allocator.Persistent);
            UnsafeUtility.MemCpy(voxs, mo.voxs, sizeof(float) * vCountXYZ);
            voxsState = (byte*)UnsafeUtility.Malloc(sizeof(byte) * vCountXYZ, UnsafeUtility.AlignOf<byte>(), Allocator.Persistent);
            UnsafeUtility.MemClear(voxsState, sizeof(byte) * vCountXYZ);

            addedRemovedVoxsI = new(mo.overlappingVoxCount, Allocator.Persistent);
            newMO = new(new(vCountXYZ), Allocator.Persistent);
            voxsToCheck = new(64, Allocator.Persistent);
            solidVoxCount = new(mo.solidVoxCount, Allocator.Persistent);
#if UNITY_EDITOR
            lToW = Matrix4x4.identity;
#endif

            int initialVerCount = mo.overlappingVoxCount * 5 * 3;
            tris = new NativeList<uint>(initialVerCount, Allocator.Persistent);
            vers = new NativeList<float3>(initialVerCount, Allocator.Persistent);
            nors = new NativeList<float3>(initialVerCount, Allocator.Persistent);
        }

        #endregion Physics
    }

    public static class MeshGenHelpers
    {
        public static T GetOrAddComponent<T>(this GameObject obj, out bool added) where T : Component
        {
            added = !obj.TryGetComponent(out T component);
            return added == false ? component : obj.AddComponent<T>();
        }

        public static T GetOrAddComponent<T>(this Component comp, out bool added) where T : Component
        {
            added = !comp.TryGetComponent(out T component);
            return added == false ? component : comp.gameObject.AddComponent<T>();
        }

        public static float Average(this Vector3 vec)
        {
            return (vec.x + vec.y + vec.z) / 3.0f;
        }

        public static float Average(this float3 vec)
        {
            return (vec.x + vec.y + vec.z) / 3.0f;
        }

#if UNITY_EDITOR
        private static readonly System.Diagnostics.Stopwatch stopwatch = new();

        public static void Debug_toggleTimer(string note = "")
        {
            if (stopwatch.IsRunning == false)
            {
                stopwatch.Restart();
            }
            else
            {
                stopwatch.Stop();
                Debug.Log(note + " time: " + stopwatch.Elapsed.TotalMilliseconds + "ms");
            }
        }
#endif
    }
}


