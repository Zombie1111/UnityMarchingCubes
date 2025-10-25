using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.Rendering;

namespace zombGen
{
    [ExecuteAlways]
    public class MarchingCubes : MonoBehaviour
    {
        [Tooltip("If null, makes the whole volume solid")]
        [SerializeField] private Mesh startShape = null;

        private void OnEnable()
        {
            Init();
        }

        private void OnDisable()
        {
            Dispose();
        }

#if UNITY_EDITOR
        private void OnValidate()
        {
            if (Application.isPlaying == true) return;
            if (lastShapeID == _currentShapeID) return;

            Dispose();
            Init();
        }
#endif

        private bool isInitilized = false;
        private VoxObject voxO;
        private int lastShapeID = 0;
        private int _currentShapeID => startShape != null ? startShape.GetInstanceID() : 0;
        private MeshCollider col = null;
        private MeshFilter mf = null;
        private MeshRenderer mr = null;
        private Mesh mesh = null;

        private void Init()
        {
            if (isInitilized == true) return;
            isInitilized = true;
            lastShapeID = _currentShapeID;
            col = this.GetOrAddComponent<MeshCollider>(out _);
            mf = this.GetOrAddComponent<MeshFilter>(out _);
            mr = this.GetOrAddComponent<MeshRenderer>(out _);
            col.sharedMesh = startShape;
            
            if (startShape != null) voxO = new(col);
            else if (TryGetComponent(out Collider othCol) == true && othCol.enabled == true)
            {
                voxO = new(col);
            }
            
            mesh = new();
            col.sharedMesh = mesh;
            mf.sharedMesh = mesh;
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            isInitilized = false;

            voxO = voxO.Dispose();
            mesh.Clear();
        }
        
        private const MeshUpdateFlags _defaultUpdateFlags =
              MeshUpdateFlags.DontRecalculateBounds
            | MeshUpdateFlags.DontValidateIndices
            | MeshUpdateFlags.DontResetBoneBounds;

        private void Update()
        {
            if (isInitilized == false) return;
            Mesh.MeshDataArray mda = Mesh.AllocateWritableMeshData(1);
            voxO.Meshify(mda[0], out Bounds boundsL);
            Mesh.ApplyAndDisposeWritableMeshData(mda, mesh, _defaultUpdateFlags);
            mesh.bounds = boundsL;
        }
    }
}


