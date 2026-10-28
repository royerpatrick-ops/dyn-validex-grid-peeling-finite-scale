"""Fourier coefficients of polygon support via analytic normal-cone integrals.

The input vertices must be in strictly convex counterclockwise order.
For coordinates z=(vertex+translation)/scale, return
  c_m=(1/(2*pi))*integral_0^(2*pi) h(phi)*exp(-1j*m*phi) dphi.
There is no angular sampling and therefore no angular-grid aliasing in
these coefficients. Evaluation uses floating-point elementary functions,
not certified interval arithmetic.

Run this file to execute finite geometry/Fourier qualification only.
"""
from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def _integrated_exponential(k,lower,width):
    """Integral exp(-i*k*phi) over [lower,lower+width], stably evaluated."""
    if k==0:return complex(width)
    return (2*math.sin(k*width/2)/k)*cmath.exp(-1j*k*(lower+width/2))


def coefficients(vertices,translation=(0.0,0.0),scale=1.0,modes=(0,1,3,6)):
    """Return {integer mode: complex coefficient}, with no angular grid.

    translation is added before division by the positive scale. The hull
    coordinates may be integers or finite floats. Adjacent repeated or
    collinear vertices, clockwise input and non-strict convexity are rejected.
    """
    if not math.isfinite(scale) or scale<=0:raise ValueError("scale must be positive and finite")
    if len(translation)!=2 or not all(math.isfinite(float(v)) for v in translation):
        raise ValueError("translation must contain two finite coordinates")
    points=[(p[0],p[1]) for p in vertices]
    if len(points)<3:raise ValueError("At least three convex vertices are required")
    if not all(math.isfinite(float(x)) and math.isfinite(float(y)) for x,y in points):
        raise ValueError("Coordinates must be finite")
    frequencies=tuple(modes)
    if any(not isinstance(m,int) for m in frequencies):raise TypeError("modes must be Python integers")
    if len(set(frequencies))!=len(frequencies):raise ValueError("modes must be distinct")
    normals=[]
    for a,b in zip(points,points[1:]+points[:1]):
        dx,dy=b[0]-a[0],b[1]-a[1]
        if dx==dy==0:raise ValueError("Repeated adjacent vertices")
        # For a CCW edge, the outward normal is the right-hand normal.
        normals.append((dy,-dx))
    widths=[];terms={m:[] for m in frequencies}
    for j,(x,y) in enumerate(points):
        previous,current=normals[j-1],normals[j]
        cross=previous[0]*current[1]-previous[1]*current[0]
        dot=previous[0]*current[0]+previous[1]*current[1]
        if cross<=0:raise ValueError("Vertices must be strictly convex and counterclockwise")
        width=math.atan2(cross,dot)
        if not 0<width<math.pi:raise ValueError("Invalid normal cone")
        widths.append(width)
        lower=math.atan2(previous[1],previous[0])
        z=complex((x+translation[0])/scale,(y+translation[1])/scale)
        for m in frequencies:
            value=(z.conjugate()*_integrated_exponential(m-1,lower,width)
                  +z*_integrated_exponential(m+1,lower,width))/(4*math.pi)
            terms[m].append(value)
    if abs(math.fsum(widths)-2*math.pi)>2e-12:
        raise ValueError("Normal cones do not make one full convex turn")
    return {m:complex(math.fsum(z.real for z in values),math.fsum(z.imag for z in values))
            for m,values in terms.items()}


def _perimeter(vertices):
    return math.fsum(math.hypot(b[0]-a[0],b[1]-a[1]) for a,b in zip(vertices,vertices[1:]+vertices[:1]))


def edge_coefficients(vertices,scale=1.0,modes=(0,3,6)):
    """Independent edge-curvature formula for modes other than +/-1.

    In distributions, h+h'' is the sum of edge lengths at their outer
    normals. Hence (1-m*m)*c_m=sum(length*exp(-i*m*normal))/(2*pi).
    These modes are unaffected by translation. This route does not use
    normal-cone integrals or their angular endpoints.
    """
    if not math.isfinite(scale) or scale<=0:raise ValueError("scale must be positive and finite")
    points=[(p[0],p[1]) for p in vertices]
    if len(points)<3:raise ValueError("At least three vertices required")
    frequencies=tuple(modes)
    if any(not isinstance(m,int) for m in frequencies):raise TypeError("modes must be Python integers")
    if any(abs(m)==1 for m in frequencies):raise ValueError("Edge formula does not determine modes +/-1")
    sums={m:[] for m in frequencies}
    for j,(x,y) in enumerate(points):
        a=points[j-1];b=points[(j+1)%len(points)]
        incoming=(x-a[0],y-a[1]);outgoing=(b[0]-x,b[1]-y)
        if incoming[0]*outgoing[1]-incoming[1]*outgoing[0]<=0:
            raise ValueError("Vertices must be strictly convex and counterclockwise")
        dx,dy=outgoing;length=math.hypot(dx,dy)
        phase=complex(dy,dx)/length
        for m in frequencies:sums[m].append((length/scale)*phase**m)
    return {m:complex(math.fsum(v.real for v in values),math.fsum(v.imag for v in values))/(2*math.pi*(1-m*m))
            for m,values in sums.items()}


def qualification():
    import numpy as np

    modes=(0,1,3,6)
    polygons={
        'square':[(-1,-1),(1,-1),(1,1),(-1,1)],
        'irregular_pentagon':[(-3,-1),(1,-2),(4,0),(2,4),(-2,3)],
    }
    checks=[]
    for n in (3,5,6,8):
        vertices=[(math.cos(2*math.pi*j/n),math.sin(2*math.pi*j/n)) for j in range(n)]
        polygons[f'regular_{n}']=vertices
        value=coefficients(vertices,modes=modes)
        for m in modes:
            if m==0:expected=n*math.sin(math.pi/n)/math.pi
            elif m%n:expected=0
            else:expected=(-1)**(m//n+1)*n*math.sin(math.pi/n)/(math.pi*(m*m-1))
            error=abs(value[m]-expected)
            assert error<2e-14,(n,m,value[m],expected)
            checks.append({'kind':'regular_polygon_closed_form','n':n,'m':m,'error':error})
    for name,vertices in polygons.items():
        original=coefficients(vertices,modes=modes)
        via_edges=edge_coefficients(vertices,modes=(0,3,6))
        for m in (0,3,6):
            error=abs(via_edges[m]-original[m])
            assert error<3e-14
            checks.append({'kind':'independent_edge_curvature_formula','polygon':name,'m':m,'error':error})
        assert abs(original[0]-_perimeter(vertices)/(2*math.pi))<3e-14
        checks.append({'kind':'perimeter_identity','polygon':name,
                       'error':abs(original[0]-_perimeter(vertices)/(2*math.pi))})
        shift=(1.25,-.375)
        shifted=coefficients(vertices,translation=shift,modes=modes)
        for m in modes:
            delta=(complex(*shift).conjugate()/2 if m==1 else 0)
            error=abs(shifted[m]-original[m]-delta)
            assert error<3e-14
            checks.append({'kind':'translation','polygon':name,'m':m,'error':error})
        scaled=coefficients(vertices,scale=2.5,modes=modes)
        for m in modes:
            error=abs(scaled[m]-original[m]/2.5)
            assert error<3e-14
            checks.append({'kind':'positive_scale','polygon':name,'m':m,'error':error})
        angle=.371
        rotated=[(math.cos(angle)*x-math.sin(angle)*y,
                  math.sin(angle)*x+math.cos(angle)*y) for x,y in vertices]
        value=coefficients(rotated,modes=modes)
        for m in modes:
            error=abs(value[m]-cmath.exp(-1j*m*angle)*original[m])
            assert error<3e-14
            checks.append({'kind':'rotation','polygon':name,'m':m,'error':error})
    quadratures=[]
    for name in ('square','regular_3','irregular_pentagon'):
        vertices=polygons[name];exact=coefficients(vertices,modes=modes)
        for n in (65536,262144):
            phi=2*math.pi*np.arange(n)/n
            arr=np.asarray(vertices)
            h=np.max(arr[:,0,None]*np.cos(phi)+arr[:,1,None]*np.sin(phi),axis=0)
            values=np.fft.rfft(h)/n
            errors={str(m):abs(complex(values[m])-exact[m]) for m in modes}
            assert max(errors.values())<2e-8
            quadratures.append({'polygon':name,'angular_nodes':n,'errors':errors})
    rejected=[]
    for vertices in ([(-1,-1),(-1,1),(1,1),(1,-1)],
                     [(-1,-1),(0,-1),(1,-1),(1,1),(-1,1)]):
        try:coefficients(vertices)
        except ValueError:rejected.append(True)
        else:raise AssertionError("Non-strict/CW polygon was accepted")
    return {'status':'PASS_ANALYTIC_POLYGON_SUPPORT_FOURIER_QUALIFICATION',
            'definition':'c_m=(1/(2*pi))*integral h(phi)*exp(-i*m*phi) dphi',
            'modes':list(modes),'angular_sampling_used_in_coefficients':False,
            'interval_certified':False,'analytic_identity_checks':len(checks),
            'max_analytic_identity_error':max(v['error'] for v in checks),
            'checks':checks,'fine_grid_cross_checks':quadratures,
            'invalid_order_or_collinear_inputs_rejected':len(rejected),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__=='__main__':
    result=qualification()
    (ROOT/'polygon_fourier_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','analytic_identity_checks',
        'max_analytic_identity_error','angular_sampling_used_in_coefficients')}))
