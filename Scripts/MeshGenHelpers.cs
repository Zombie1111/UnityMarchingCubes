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
    public readonly unsafe struct MarchingObject
    {
        public readonly float* voxs;
        public readonly float3 start;
        public readonly float3 voxSize;
        public readonly float surface;

        public readonly NativeList<uint> tris;
        public readonly NativeList<float3> vers;
        public readonly NativeList<float3> nors;

        /// <summary>
        /// Lenght of voxs, -1 if invalid
        /// </summary>
        public readonly int vCountXYZ;
        public readonly int vCountYZ;
        public readonly int vCountZ;
        public readonly int vCountY;
        public readonly int vCountX;

        #region Voxelization

        /// <summary>
        /// DO NOT ACCESS this after disposing (Returns new cleared VoxObject)
        /// </summary>
        public MarchingObject Dispose()
        {
            if (vCountXYZ < 0) return this;

            UnsafeUtility.Free(voxs, Allocator.Persistent);
            tris.Dispose();
            vers.Dispose();
            nors.Dispose();
            return new(true);
        }

        private MarchingObject(bool unused)
        {
            voxs = null;
            start = Vector3.zero;
            voxSize = Vector3.zero;
            surface = 0.0f;

            tris = new();
            vers = new();
            nors = new();

            vCountXYZ = -1;
            vCountYZ = -1;
            vCountZ = -1;
            vCountY = -1;
            vCountX = -1;
        }

        /// <summary>
        /// Voxelizes the given collider. (Mainthread only)
        /// DONT FORGET TO CALL this.Dispose() in OnDestroy or similar
        /// </summary>
        public MarchingObject(Collider col, bool tryFillConvaveInteriors = true, bool tryStripOuterLayer = false,
            float surface = 0.0f, float isoLevel = 1.0f)
        {
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

            Vector3 bStart = colBounds.min - voxelSize;
            Vector3 bMax = colBounds.max;
            int ogLayer = col.gameObject.layer;
            col.gameObject.layer = MeshGenGlobals.voxelizeTempLayer;
            if (ogLayer == MeshGenGlobals.voxelizeTempLayer) Debug.LogWarning(col.transform.name + " already had voxelizeTempLayer layer, not recommended");
            LayerMask layerMask = 1 << MeshGenGlobals.voxelizeTempLayer;
            PhysicsScene colPhyScene = col.gameObject.scene.GetPhysicsScene();

            //Overlap tests
            Vector3 voxelSizeHalf = voxelSize * 0.501f;
            Vector3 voxelSizeHalfReal = voxelSize * 0.5f;
            int vCountXYZ = vCountX * vCountY * vCountZ;
            float* voxs = null;
            int overlappingVoxCount = 0;

            if (vCountXYZ <= 0)
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

            voxs = (float*)UnsafeUtility.Malloc(vCountXYZ * sizeof(float),
                UnsafeUtility.AlignOf<float>(), Allocator.Persistent);

            Parallel.For(0, vCountXYZ, vI =>
            {
                int maxRI = (vI * MeshGenGlobals.voxelizeMaxHits) + MeshGenGlobals.voxelizeMaxHits;
                for (int rI = vI * MeshGenGlobals.voxelizeMaxHits; rI < maxRI; rI++)
                {
                    if (results[rI].instanceID != colId) continue;
                    if (vI < 0 || vI >= vCountXYZ)
                    {
                        continue;
                    }

                    voxs[vI] = isoLevel;
                    overlappingVoxCount++;
                    return;
                }

                voxs[vI] = -1.0f;
            });

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
                    if (voxs[vI] == -1.0f) voxs[vI] = isoLevel;
                    else if (voxs[vI] == 3.0f) voxs[vI] = -1.0f;
                }
            }

        SkipCreation:;

            this.voxs = voxs;
            start = colWToL.MultiplyPoint3x4(bStart);
            voxSize = inverseScale * MeshGenGlobals.voxelSizeWorld;
            this.vCountXYZ = vCountXYZ;
            this.vCountYZ = vCountYZ;
            this.vCountZ = vCountZ;
            this.vCountY = vCountY;
            this.vCountX = vCountX;
            this.surface = surface;

            int initialVerCount = overlappingVoxCount * 5 * 3;
            tris = new NativeList<uint>(initialVerCount, Allocator.Persistent);
            vers = new NativeList<float3>(initialVerCount, Allocator.Persistent);
            nors = new NativeList<float3>(initialVerCount, Allocator.Persistent);

            col.gameObject.layer = ogLayer;
            col.transform.rotation = ogRot;
        }
        #endregion Voxelization

        #region Meshification

        public void Meshify(Mesh.MeshData md, out Bounds boundsLocal)
        {
            //Alloc
            float* cubeValues = (float*)UnsafeUtility.Malloc(8 * sizeof(float),
                UnsafeUtility.AlignOf<float>(), Allocator.Temp);
            float3* edgeVertex = (float3*)UnsafeUtility.Malloc(12 * 3 * sizeof(float),
                UnsafeUtility.AlignOf<float3>(), Allocator.Temp);

            //We could persistent allocate these in a readonly struct and pass along that, worth it?
            int* offsets = MeshGenLookup.GetOffsets(vCountZ, vCountYZ);
            int* cubeEdgeFlags = MeshGenLookup.GetCubeEdgeFlags();
            int2* edgeConnections = MeshGenLookup.GetEdgeConnections();
            float3* edgeDirections = MeshGenLookup.GetEdgeDirections();
            float3* vertexOffsets = MeshGenLookup.GetVertexOffsets();
            int* triangleConnectionTables = MeshGenLookup.GetTriangleConnectionTables();

            float3 voxSize = this.voxSize;
            float3 minPos = new(69420.0f, 69420.0f, 69420.0f);
            float3 maxPos = -minPos;
            vers.Clear();
            tris.Clear();
            nors.Clear();

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
                    //int remainderAfterZ = nVI % vCountYZ;
                    //int x = nVI / vCountYZ;
                    //int y = remainderAfterZ / vCountZ;
                    //int z = remainderAfterZ % vCountZ;
                    //if (x < 0 || x > vMaxX
                    //    || y < 0 || y > vMaxY
                    //    || z < 0 || z > vMaxZ)
                    if (IsVoxIndexInBounds(nVI) == false)
                    {//Bounds check
                        cubeValues[i] = -1.0f;
                        continue;
                    }

                    cubeValues[i] = voxs[nVI];
                    if (voxs[nVI] >= surface) cubeIndex |= 1 << i;
                }

                int edgeFlag = cubeEdgeFlags[cubeIndex];
                if (edgeFlag == 0) continue;

                for (int i = 0; i < 12; i++)
                {
                    //if there is an intersection on this edge
                    if ((edgeFlag & (1 << i)) == 0) continue;
                    int2 edgeCon = edgeConnections[i];
                    float offset = GetOffset(cubeValues[edgeCon.x],
                        cubeValues[edgeCon.y]);

                    edgeVertex[i] = posL + ((vertexOffsets[edgeCon.x] + offset * edgeDirections[i]) * voxSize);
                }

                int cubeIndex_tri = cubeIndex * 16;
                for (int i = 0; i < 16; i += 3)
                {
                    if (triangleConnectionTables[cubeIndex_tri + i] < 0) break;

                    float3 ev0 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i]];
                    float3 ev1 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i + 1]];
                    float3 ev2 = edgeVertex[triangleConnectionTables[cubeIndex_tri + i + 2]];
                    minPos = math.min(minPos, math.min(ev0, math.min(ev1, ev2)));
                    maxPos = math.max(maxPos, math.max(ev0, math.max(ev1, ev2)));
                    vers.Add(ev0);
                    vers.Add(ev1);
                    vers.Add(ev2);

                    Vector3 nor = math.normalize(math.cross(ev1 - ev2, ev0 - ev2));
                    nors.Add(nor);//Cross order is reversed as tri order
                    nors.Add(nor);
                    nors.Add(nor);

                    uint idx = (uint)tris.Length;
                    tris.Add(idx + 2);
                    tris.Add(idx + 1);
                    tris.Add(idx);
                }
            }

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
            UnsafeUtility.Free(cubeEdgeFlags, Allocator.Temp);
            UnsafeUtility.Free(edgeConnections, Allocator.Temp);
            UnsafeUtility.Free(edgeDirections, Allocator.Temp);
            UnsafeUtility.Free(vertexOffsets, Allocator.Temp);
            UnsafeUtility.Free(triangleConnectionTables, Allocator.Temp);
        }

        #endregion Meshification

        #region Base Helpers

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
                Mathf.Clamp(posL.x, voxSize.x * 1.5f * voxMargin, (vCountX * voxSize.x) - (voxSize.x * 1.5f * voxMargin)),
                Mathf.Clamp(posL.y, voxSize.y * 1.5f * voxMargin, (vCountY * voxSize.y) - (voxSize.y * 1.5f * voxMargin)),
                Mathf.Clamp(posL.z, voxSize.z * 1.5f * voxMargin, (vCountZ * voxSize.z) - (voxSize.z * 1.5f * voxMargin))
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

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private readonly float GetOffset(float v1, float v2)
        {
            float delta = v2 - v1;
            return (delta == 0.0f) ? surface : (surface - v1) / delta;
        }

        #endregion Base Helpers

        #region Methods

        public readonly void SetVoxelsAround(int voxI, float radiusLocal, float newValue = -1.0f)
        {
            int maxX = Mathf.CeilToInt(radiusLocal * voxSize.x);
            int maxY = Mathf.CeilToInt(radiusLocal * voxSize.y);
            int maxZ = Mathf.CeilToInt(radiusLocal * voxSize.z);
            float radiusSQ = radiusLocal * radiusLocal;
            int margin = newValue >= surface ? 1 : 0;
            //Margin when creating as marching cubes requires 1 row of empty voxels to work properly

            for (int z = -maxZ; z <= maxZ; z++)
            {
                for (int y = -maxY; y <= maxY; y++)
                {
                    for (int x = -maxX; x <= maxX; x++)
                    {
                        float dxw = x * voxSize.x;
                        float dyw = y * voxSize.y;
                        float dzw = z * voxSize.z;
                        if (dxw * dxw + dyw * dyw + dzw * dzw > radiusSQ) continue;

                        int vI = voxI + z + (vCountZ * y) + (vCountYZ * x);
                        if (IsVoxIndexInBounds(vI, margin) == false) continue;
                        voxs[vI] = newValue;
                    }
                }
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly void RemoveVoxelsAround(float3 posL, float radiusLocal)
        {
            SetVoxelsAround(PosLToVoxIndex(posL), radiusLocal, -1.0f);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public readonly void CreateVoxelsAround(float3 posL, float radiusLocal)
        {
            SetVoxelsAround(PosLToVoxIndex(posL), radiusLocal, 1.0f);
        }
        #endregion Methods
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
    }
}


