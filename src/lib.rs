pub use kurbo;
use kurbo::{Point, Rect, Vec2, PathEl, BezPath, ParamCurve};

#[derive(Debug, Copy, Clone)]
pub enum CubicApprox {
    //
    Linear,
    // Generates line segments
    Flatten(f64),
    // Fast but rough approximation (1-2 quads)
    Midpoint,
    // Generates quadratic segments
    Lyon(f64),
}

pub fn break_path(orig: &BezPath, cubic_approx: CubicApprox) -> BezPath {
    let mut initial = Point::ORIGIN;
    let mut p0 = Point::ORIGIN;
    let mut path = BezPath::new();

    fn add_quad(path: &mut BezPath, p0: Point, p1: Point, p2: Point) {
        // Split quadratic at parametric point t
        //
        // Returns new control points and midpoint
        fn split(t: f64, p0: Point, p1: Point, p2: Point) -> [Point; 3] {
            let pa = p0.lerp(p1, t);
            let pc = p1.lerp(p2, t);
            let pb = pa.lerp(pc, t);

            [pa, pb, pc]
        }

        let aabb = Rect::from_points(p0, p1);
        let tx = if aabb.min_x() > p1.x || p1.x > aabb.max_x() {
            Some((p0.x - p1.x) / (p0.x - 2.0 * p1.x + p2.x))
        } else {
            None
        };
        let ty = if aabb.min_y() > p1.y || p1.y > aabb.max_y() {
            Some((p0.y - p1.y) / (p0.y - 2.0 * p1.y + p2.y))
        } else {
            None
        };

        match (tx, ty) {
            (Some(tx), Some(ty)) => {
                let t0 = tx.min(ty);
                let t1 = (tx.max(ty) - t0) / (1.0 - t0);

                let [pa0, pb0, pc0] = split(t0, p0, p1, p2);
                let [pa1, pb1, pc1] = split(t1, pb0, pc0, p2);

                path.quad_to(pa0, pb0);
                path.quad_to(pa1, pb1);
                path.quad_to(pc1, p2);
            }
            (Some(t), None) | (None, Some(t)) => {
                let [pa, pb, pc] = split(t, p0, p1, p2);
                path.quad_to(pa, pb);
                path.quad_to(pc, p2);
            }
            (None, None) => {
                path.quad_to(p1, p2);
            }
        }
    }

    for elem in orig {
        match elem {
            PathEl::MoveTo(p) => {
                path.move_to(p);
                p0 = p;
                initial = p;
            }
            PathEl::LineTo(p) => {
                path.line_to(p);
                p0 = p;
            }
            PathEl::QuadTo(p1, p2) => {
                add_quad(&mut path, p0, p1, p2);
                p0 = p2;
            }
            PathEl::CurveTo(p1, p2, p3) => {
                match cubic_approx {
                    CubicApprox::Linear => {
                        path.line_to(p3);
                    }
                    CubicApprox::Flatten(tolerance) => {
                        let mut subpath = BezPath::new();
                        subpath.move_to(p0);
                        subpath.curve_to(p1, p2, p3);
                        kurbo::flatten(subpath, tolerance, |el| match el {
                            PathEl::MoveTo(_) => {}
                            PathEl::LineTo(p) => {
                                path.line_to(p);
                                p0 = p;
                            }
                            _ => unreachable!(),
                        });
                    }
                    CubicApprox::Midpoint => {
                        // 3.5 Alternative approximation of cubic curves
                        if p0 == p1 {
                            add_quad(&mut path, p0, p2, p3);
                        } else if p2 == p3 {
                            add_quad(&mut path, p0, p1, p3);
                        } else {
                            let p_ca = p0.lerp(p1, 0.75);
                            let p_cb = p3.lerp(p2, 0.75);
                            let p_m = p_ca.midpoint(p_cb);
                            add_quad(&mut path, p0, p_ca, p_m);
                            add_quad(&mut path, p_m, p_cb, p3);
                        }
                    }
                    CubicApprox::Lyon(tolerance) => {
                        use lyon_geom::{
                            cubic_bezier::CubicBezierSegment,
                            cubic_to_quadratic::cubic_to_quadratics,
                        };

                        // monotonic variant appears to be buggy (v0.15)
                        cubic_to_quadratics(
                            &CubicBezierSegment {
                                from: [p0.x, p0.y].into(),
                                ctrl1: [p1.x, p1.y].into(),
                                ctrl2: [p2.x, p2.y].into(),
                                to: [p3.x, p3.y].into(),
                            },
                            tolerance,
                            &mut |segment| {
                                add_quad(
                                    &mut path,
                                    Point::new(segment.from.x, segment.from.y),
                                    Point::new(segment.ctrl.x, segment.ctrl.y),
                                    Point::new(segment.to.x, segment.to.y),
                                );
                            },
                        );
                    }
                }

                p0 = p3;
            }
            PathEl::ClosePath => {
                path.close_path();
                p0 = initial;
            }
        }
    }
    path
}

pub fn monotonize_quads(orig: &BezPath) -> BezPath {
    let mut initial = Point::ORIGIN;
    let mut p0 = Point::ORIGIN;
    let mut path = BezPath::new();

    fn split_quad(path: &mut BezPath, p0: Point, p1: Point, p2: Point) {
        let min = Point::new(p0.x.min(p2.x), p0.y.min(p2.y));
        let max = Point::new(p0.x.max(p2.x), p0.y.max(p2.y));

        let tx = if p1.x < min.x || max.x < p1.x {
            Some((p0.x - p1.x) / (p0.x - 2.0 * p1.x + p2.x))
        } else {
            None
        };

        let ty = if p1.y < min.y || max.y < p1.y {
            Some((p0.y - p1.y) / (p0.y - 2.0 * p1.y + p2.y))
        } else {
            None
        };

        match (tx, ty) {
            (Some(tx), Some(ty)) => {
                let (t0, t1) = (tx.min(ty), tx.max(ty));

                let t = t0;
                let c2 = kurbo::QuadBez::new(p0, p1, p2).eval(t);
                let c1 = kurbo::Line::new(p0, p1).eval(t);
                let c3 = kurbo::Line::new(p1, p2).eval(t);

                let t = (t1 - t0) / (1.0 - t1);
                let c5 = kurbo::QuadBez::new(c2, c3, p2).eval(t);
                let c4 = kurbo::Line::new(c2, c3).eval(t);
                let c6 = kurbo::Line::new(c3, p2).eval(t);

                path.quad_to(c1, c2);
                path.quad_to(c4, c5);
                path.quad_to(c6, p2);
            }
            (Some(t), None) | (None, Some(t)) => {
                let p = kurbo::QuadBez::new(p0, p1, p2).eval(t);
                let p10 = kurbo::Line::new(p0, p1).eval(t);
                let p11 = kurbo::Line::new(p1, p2).eval(t);
                path.quad_to(p10, p);
                path.quad_to(p11, p2);
            }
            (None, None) => path.quad_to(p1, p2),
        }
    }

    for elem in orig {
        match elem {
            PathEl::MoveTo(p) => {
                path.move_to(p);
                p0 = p;
                initial = p;
            }
            PathEl::LineTo(p) => {
                path.line_to(p);
                p0 = p;
            }
            PathEl::QuadTo(p1, p2) => {
                split_quad(&mut path, p0, p1, p2);
                p0 = p2;
            }
            PathEl::CurveTo(p1, p2, p3) => {
                // quads only
                path.curve_to(p1, p2, p3);
                p0 = p3;
            }
            PathEl::ClosePath => {
                path.close_path();
                p0 = initial;
            }
        }
    }

    path
}
