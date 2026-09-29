# Cylinder geometry (user PhD mesh)

Canonical file: `cylinder.msh`.

Physical patches after import:
- **top** — CAP/PAM exposure (`fixedValue` H₂O₂ / NO₂⁻)
- **bottom** — case-specific scalar boundary condition
- **walls** — case-specific scalar boundary condition

Domain (from mesh): radius 0.5 mm, height 1 mm.
Bounding box: `(-0.0005,-0.0005,0)` → `(0.0005,0.0005,0.001)`.

Re-import into a case:
```bash
./import_cylinder_mesh.sh ../../participants/tissue_41_diffusion
```

The importer converts through a temporary Gmsh v2 file and deletes it.
`constant/cylinder.msh` therefore remains the only `.msh` file in each case.
