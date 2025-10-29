//Copyright 2025 David Westberg (MIT) https://github.com/Zombie1111/UnityMarchingCubes
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

        //Maximum volume of a chunk for it to disconnect is maxChunkCheckRadius^3
        public const int maxChunkCheckRadius = 20;

        public const byte kinVoxFlag = 1 << 7;
        public const byte delVoxFlag = 1 << 6;
        public const byte kinDelVoxFlag = delVoxFlag | kinVoxFlag;
    }
}

