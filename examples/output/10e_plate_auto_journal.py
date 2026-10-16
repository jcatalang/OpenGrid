import opengrid

og = opengrid.Session()
og.create_line(p0=(0.0, 0.0, 0.0), p1=(4.0, 0.0, 0.0), name='b0')
og.create_line(p0=(4.0, 0.0, 0.0), p1=(4.0, 2.0, 0.0), name='b1')
og.create_line(p0=(4.0, 2.0, 0.0), p1=(0.0, 2.0, 0.0), name='b2')
og.create_line(p0=(0.0, 2.0, 0.0), p1=(0.0, 0.0, 0.0), name='b3')
og.create_arc(center=(1.0, 1.0, 0.0), radius=0.3, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), start_angle=0.0, end_angle=6.283185307179586, name='hole0')
og.create_arc(center=(2.0, 0.8, 0.0), radius=0.3, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), start_angle=0.0, end_angle=6.283185307179586, name='hole1')
og.create_arc(center=(3.0, 1.2, 0.0), radius=0.3, normal=(0.0, 0.0, 1.0), x_axis=(1.0, 0.0, 0.0), start_angle=0.0, end_angle=6.283185307179586, name='hole2')
og.auto_topology(curves=('b0', 'b1', 'b2', 'b3', 'hole0', 'hole1', 'hole2'))
