using System;
using System.Collections;
using System.Collections.Generic;
using System.Threading.Tasks;
using Unity.Collections;
using UnityEngine;

namespace zombGen
{
    public readonly struct VoxObject
    {
        public readonly byte[] voxs;
        public readonly Vector3 start;
        public readonly Vector3 voxSize;
        public readonly float surface;

        public readonly int vCountXYZ;
        public readonly int vCountYZ;
        public readonly int vCountZ;
        public readonly int vCountY;
        public readonly int vCountX;

        #region Voxelization

        /// <summary>
        /// Voxelizes the given collider. The voxels has the size defined in voxGlobalSettings.cs
        /// </summary>
        public VoxObject(Collider col, bool tryFillConvaveInteriors = true)
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
            byte[] voxs = new byte[vCountXYZ];
            if (vCountXYZ <= 0)
            {
                Debug.LogWarning(col.transform.name + " wont contain any voxels! Is the object too small?");
                goto SkipCreation;
            }

            if (vCountXYZ > MeshGenGlobals.maxVoxelsInExtent
                * MeshGenGlobals.maxVoxelsInExtent
                * MeshGenGlobals.maxVoxelsInExtent)
            {
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

            Parallel.For(0, vCountXYZ, vI =>
            {
                int maxRI = (vI * MeshGenGlobals.voxelizeMaxHits) + MeshGenGlobals.voxelizeMaxHits;
                for (int rI = vI * MeshGenGlobals.voxelizeMaxHits; rI < maxRI; rI++)
                {
                    if (results[rI].instanceID != colId) continue;
                    voxs[vI] = 1;
                    return;
                }
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
                    if (voxs[vI] != 0) continue;
                    voxs[vI] = 2;
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
                        if (voxs[vI] != 2) continue;
                        voxs[vI] = 3;

                        for (int i = 0; i < 6; i++)
                        {
                            int nVI = vI + offsets[i];
                            int remainderAfterZ = nVI % vCountYZ;
                            int x = nVI / vCountYZ;
                            int y = remainderAfterZ / vCountZ;
                            int z = remainderAfterZ % vCountZ;
                            if ((x < 0 || x > vMaxX ? 1 : 0)//Bounds check
                            + (y < 0 || y > vMaxY ? 1 : 0)
                            + (z < 0 || z > vMaxZ ? 1 : 0) > 0) continue;

                            if (voxs[nVI] != 0) continue;
                            voxs[nVI] = 2;
                            spreadedAny = true;
                        }
                    }
                }

                for (int vI = 0; vI < vCountXYZ; vI++)
                {
                    if (voxs[vI] == 0) voxs[vI] = 1;
                    else if (voxs[vI] == 3) voxs[vI] = 0;
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
            surface = 0.0f;

            col.gameObject.layer = ogLayer;
            col.transform.rotation = ogRot;
        }
        #endregion Voxelization

        public void Meshify(Transform ttt)
        {
            Vector3 scale = ttt.lossyScale;
            Vector3 voxSize = new(this.voxSize.x, this.voxSize.y, this.voxSize.z);
            List<Vector3> vers = new(64);
            List<int> tris = new(64);
            float[] cubeValues = new float[8];

            var offsets = MeshGenLookup.GetOffsets(vCountZ, vCountYZ);
            var edgeVertex = new Vector3[12];

            for (int vI = 0; vI < vCountXYZ; vI++)
            {
                Vector3 posL = VoxIndexToPosL(vI);
                int cubeIndex = 0;
                int vMaxX = vCountX - 1;
                int vMaxY = vCountY - 1;
                int vMaxZ = vCountZ - 1;

                for (int i = 0; i < 8; i++)
                {
                    int nVI = vI + offsets[i];
                    int remainderAfterZ = nVI % vCountYZ;
                    int x = nVI / vCountYZ;
                    int y = remainderAfterZ / vCountZ;
                    int z = remainderAfterZ % vCountZ;
                    if ((x < 0 || x > vMaxX ? 1 : 0)//Bounds check
                    + (y < 0 || y > vMaxY ? 1 : 0)
                    + (z < 0 || z > vMaxZ ? 1 : 0) > 0
                    || voxs[nVI] == 0)
                    {
                        cubeValues[i] = -1.0f;
                        continue;
                    }

                    cubeValues[i] = 1.0f;
                    cubeIndex |= 1 << i;
                }

                int edgeFlag = MeshGenLookup.CubeEdgeFlags[cubeIndex];
                if (edgeFlag == 0) continue;

                for (int i = 0; i < 12; i++)
                {
                    //if there is an intersection on this edge
                    if ((edgeFlag & (1 << i)) != 0)
                    {
                        float offset = GetOffset(cubeValues[MeshGenLookup.EdgeConnection[i, 0]],
                            cubeValues[MeshGenLookup.EdgeConnection[i, 1]]);

                        edgeVertex[i].x = posL.x
                            + ((MeshGenLookup.VertexOffset[MeshGenLookup.EdgeConnection[i, 0], 0]
                            + offset * MeshGenLookup.EdgeDirection[i, 0]) * voxSize.x);
                        edgeVertex[i].y = posL.y
                            + ((MeshGenLookup.VertexOffset[MeshGenLookup.EdgeConnection[i, 0], 1]
                            + offset * MeshGenLookup.EdgeDirection[i, 1]) * voxSize.y);
                        edgeVertex[i].z = posL.z
                            + ((MeshGenLookup.VertexOffset[MeshGenLookup.EdgeConnection[i, 0], 2]
                            + offset * MeshGenLookup.EdgeDirection[i, 2]) * voxSize.z);
                    }
                }

                for (int i = 0; i < 5; i++)
                {
                    if (MeshGenLookup.TriangleConnectionTable[cubeIndex, 3 * i] < 0) break;
                    int idx = vers.Count;

                    for (int ii = 0; ii < 3; ii++)
                    {
                        int verI = MeshGenLookup.TriangleConnectionTable[cubeIndex, 3 * i + ii];
                        tris.Add(idx + (2 - ii));//Invert ii order to flip faces
                        vers.Add(edgeVertex[verI]);
                    }
                }
            }

            Mesh m = new()
            {
                vertices = vers.ToArray(),
                triangles = tris.ToArray(),
            };

            m.RecalculateNormals();
            m.RecalculateBounds();
            Gizmos.DrawMesh(m, 0, ttt.position, ttt.rotation, ttt.lossyScale);
        }

        public readonly Vector3 VoxIndexToPosL(int voxI)
        {
            int remainderAfterZ = voxI % vCountYZ;
            return new Vector3((voxI / vCountYZ) * voxSize.x,
                (remainderAfterZ / vCountZ) * voxSize.y,
                (remainderAfterZ % vCountZ) * voxSize.z) + start;
        }

        public readonly float GetOffset(float v1, float v2)
        {
            float delta = v2 - v1;
            return (delta == 0.0f) ? surface : (surface - v1) / delta;
        }
    }

    public static class MeshGenHelpers
    {
        #region Generic

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

        #endregion Generic

       
    }
}


