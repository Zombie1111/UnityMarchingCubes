using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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

        private void OnDrawGizmosSelected()
        {
            if (voxO.vCountXYZ == 0)
            {
                if (startShape != null) Gizmos.DrawWireMesh(startShape, 0, transform.position, transform.rotation, transform.localScale);
                return;
            }

            voxO.Meshify(transform);
        }

        private VoxObject voxO;
        private int lastShapeID = 0;
        private int _currentShapeID => startShape != null ? startShape.GetInstanceID() : 0;
        private MeshCollider col = null;

        private void Init()
        {
            if (isInitilized == true) return;
            isInitilized = true;
            lastShapeID = _currentShapeID;
            col = this.GetOrAddComponent<MeshCollider>(out _);

            if (startShape != null)
            {
                col.sharedMesh = startShape;
                voxO = new(col);
            }
            else if (TryGetComponent(out Collider othCol) == true && othCol.enabled == true)
            {
                voxO = new(col);
            }
        }

        private void Dispose()
        {
            if (isInitilized == false) return;
            isInitilized = false;
        }
    }
}


