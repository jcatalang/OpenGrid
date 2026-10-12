import opengrid

og = opengrid.Session()
og.create_arc(center=(0.0, 0.0, 0.0), radius=0.5, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), start_angle=0.0, end_angle=6.283185307179586, name='circle0')
og.create_arc(center=(0.0, 0.0, 1.0), radius=0.5, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), start_angle=0.0, end_angle=6.283185307179586, name='circle1')
og.create_line(p0=(0.5, 0.0, 0.0), p1=(0.5, 0.0, 1.0), name='seam')
og.create_connector(segments=({'curve': 'circle0', 't0': 0.0, 't1': 1.0},), dimension=129, mode='arc_length', distribution={'type': 'uniform', 'begin_spacing': None, 'end_spacing': None, 'growth_ratio': None}, breakpoints=(), propagate_spacing=True, name='c-circle0')
og.create_connector(segments=({'curve': 'circle1', 't0': 0.0, 't1': 1.0},), dimension=129, mode='arc_length', distribution={'type': 'uniform', 'begin_spacing': None, 'end_spacing': None, 'growth_ratio': None}, breakpoints=(), propagate_spacing=True, name='c-circle1')
og.create_connector(segments=({'curve': 'seam', 't0': 0.0, 't1': 1.0},), dimension=21, mode='arc_length', distribution={'type': 'uniform', 'begin_spacing': None, 'end_spacing': None, 'growth_ratio': None}, breakpoints=(), propagate_spacing=True, name='c-seam')
og.create_domain(edges=(({'connector': 'c-circle0', 'reversed': False},), ({'connector': 'c-seam', 'reversed': False},), ({'connector': 'c-circle1', 'reversed': True},), ({'connector': 'c-seam', 'reversed': True},)), generator={'method': 'tfi', 'hermite_walls': (0,), 'hermite_strength': 1.0, 'elliptic': None}, surface=None, name='wall')
og.create_block(faces=None, extrusion={'domain': 'wall', 'mode': 'hyperbolic', 'vector': (0.0, 0.0, 1.0), 'axis_point': (0.0, 0.0, 0.0), 'angle': 0.0, 'distance': 10.0, 'dimension': 61, 'distribution': {'type': 'tanh', 'begin_spacing': 0.0001, 'end_spacing': None, 'growth_ratio': None}, 'edges': ('periodic', 'periodic', 'constant_z', 'constant_z'), 'dissipation': 0.2}, smoothing=None, quality_floor={'min_jacobian': 0.05, 'max_skewness': 0.95}, name='ogrid')
og.create_patch(patch_type='wall', members=('wall',), periodic=None, name='cylinder')
og.create_patch(patch_type='farfield', members=({'block': 'ogrid', 'face': 'kmax'},), periodic=None, name='farfield')
og.create_patch(patch_type='symmetry', members=({'block': 'ogrid', 'face': 'jmin'}, {'block': 'ogrid', 'face': 'jmax'}), periodic=None, name='ends')
