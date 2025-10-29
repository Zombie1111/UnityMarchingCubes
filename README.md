<h1>Unity Marching Cubes</h1>

## Overview

Burst and job compatible implementation of Marching Cubes

<img src="https://media4.giphy.com/media/v1.Y2lkPTc5MGI3NjExN3h0MjIzMXhvMmU4bWoydTdleDY0MW01bndnc2F3aGdwMWo4YXg1eiZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/6qtgyiIbrfH5IwUrob/giphy.gif" width="100%" height="100%"/>

## Key Features
<ul>
<li>Multithreaded (Burst and Job)</li>
<li>Smooth normals</li>
<li>Create a marched object from any mesh</li>
<li>Collision (Mesh collider baked in job or many primitive shapes)</li>
<li>Detects loose chunks and splits them into seperate objects with rigidbody</li>
</ul>

## Instructions
**Requirements** (Should work in other versions and pipelines)
<ul>
<li>Unity 2023.2.20f1 (Built-in)</li>
<li>Burst 1.8.13</li>
<li>Collections 2.1.4</li>
<li>Mathematics 1.2.6</li>
</ul>

**General Setup**
<ol>
  <li>Download and copy the <code>Demo</code> and <code>Scripts</code> folders into an empty folder inside your <code>Assets</code> directory</li>
  <li>Add <code>MeshGenObject.cs</code> script to any empty Gameobject and assign a mesh (preferable watertight) to the <code>StartShape</code> field</li>
  <li>See <code>Demo/meshGenDemo</code> scene for example</li>
</ol>

## License
MIT - See the `LICENSE` file for more details.