using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace zombGen
{
    public static class MeshGenGlobals
    {
        public const float voxelSizeWorld = 0.5f;

        //Unity crashes if more overlaps that this occures!!!!
        public const int voxelizeMaxHits = 4;

        //Recommended that no other colliders are on this layer, to prevent crash explained above
        public const int voxelizeTempLayer = 2;

        //Maximum allowed size of object is maxVoxelsInExtent * voxelSizeWorld
        public const int maxVoxelsInExtent = 128;
    }
}

