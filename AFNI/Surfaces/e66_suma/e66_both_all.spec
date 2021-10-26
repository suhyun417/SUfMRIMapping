# define the group
  Group = e66

# define various States
  StateDef = smoothwm
  StateDef = inflated_lh
  StateDef = inflated_rh
  StateDef = veryinflated_lh
  StateDef = veryinflated_rh
  StateDef = flat_lh
  StateDef = flat_rh

NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_LH_FS_surf.asc
  LocalDomainParent = SAME
  SurfaceState = smoothwm
  EmbedDimension = 3
  Anatomical = Y

NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_RH_FS_surf.asc
  LocalDomainParent = SAME
  SurfaceState = smoothwm
  EmbedDimension = 3
  Anatomical = Y

NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_LH_FS_inflated.asc
  LocalDomainParent = e66_LH_FS_surf.asc
  SurfaceState = inflated_lh
  EmbedDimension = 3
  Anatomical = N


NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_RH_FS_inflated.asc
  LocalDomainParent = e66_RH_FS_surf.asc
  SurfaceState = inflated_rh
  EmbedDimension = 3
  Anatomical = N

NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_LH_FS_veryinflated.asc
  LocalDomainParent = e66_LH_FS_surf.asc
  SurfaceState = veryinflated_lh
  EmbedDimension = 3
  Anatomical = N


NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_RH_FS_veryinflated.asc
  LocalDomainParent = e66_RH_FS_surf.asc
  SurfaceState = veryinflated_rh
  EmbedDimension = 3
  Anatomical = N

NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_LH_FS_flat.asc
  LocalDomainParent = e66_LH_FS_surf.asc
  SurfaceState = flat_lh
  EmbedDimension = 3
  Anatomical = N


NewSurface
  SurfaceFormat = ASCII
  SurfaceType = FreeSurfer
  FreeSurferSurface = e66_RH_FS_flat.asc
  LocalDomainParent = e66_RH_FS_surf.asc
  SurfaceState = flat_rh
  EmbedDimension = 3
  Anatomical = N
