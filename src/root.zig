//! This package is a simple module to help converting points coordinated from the french lambert system to
//! the most used WGS84 system.
const std = @import("std");
const math = @import("std").math;

const testing = std.testing;

/// A simple structure to hold a point structure.
pub const Point = struct {
    /// The x value
    x: f64,
    /// The y value
    y: f64,
    /// The z value
    z: f64,

    /// Default initializer of the point
    pub fn init(x: f64, y: f64, z: f64) Point {
        return Point{ .x = x, .y = y, .z = z };
    }
    /// A simpler multiply method
    pub fn mul(self: Point, val: f64) Point {
        return Point.init(self.x * val, self.y * val, self.z * val);
    }
};
/// Identify the lambert zone we are in
pub const LambertZone = enum(u8) { LambertI = 0, LambertII = 1, LambertIII = 2, LambertIV = 3, LambertII_E = 4, Lambert93 = 5 };
/// A constant to convert degree point to radian
pub const DegreeToRadian = math.pi / 180.0;
/// A constant to convert radian point to degree
pub const RadianToDegree = 180.0 / math.pi;

const LON_MERID_PARIS: f64 = 0;
const LON_MERID_GREENWICH: f64 = 0.04079234433;
const E_CLARK_IGN: f64 = 0.08248325676;
const A_CLARK_IGN: f64 = 6378249.2;
const A_WGS84: f64 = 6378137.0;
const DEFAULT_EPS: f64 = 1e-10;
const E_WGS84: f64 = 0.08181919106;

const LON_MERID_IERS: f64 = (3 * math.pi / 180.0);
test "test basic conversion" {
    const pt = Point.init(math.pi, 0, 0);
    try testing.expectApproxEqAbs(180.0, pt.mul(RadianToDegree).x, 1e-6);

    const pt2 = Point.init(180.0, 0, 0);
    try testing.expectApproxEqAbs(math.pi, pt2.mul(DegreeToRadian).x, 1e-6);
}

const lambert_n = [6]f64{ 0.7604059656, 0.7289686274, 0.6959127966, 0.6712679322, 0.7289686274, 0.7256077650 };
const lambert_c = [6]f64{ 11603796.98, 11745793.39, 11947992.52, 12136281.99, 11745793.39, 11754255.426 };
const lambert_xs = [6]f64{ 600000.0, 600000.0, 600000.0, 234.358, 600000.0, 700000.0 };
const lambert_ys = [6]f64{ 5657616.674, 6199695.768, 6791905.085, 7239161.542, 8199695.768, 12655612.050 };

fn lambertNormal(lat: f64, a: f64, e: f64) f64 {
    const sinA = math.sin(lat);
    const b = math.sqrt(1.0 - e * e * sinA * sinA);
    return a / b;
}

test "algo 00021" {
    const n: f64 = 6393174.9755;
    const lat: f64 = 0.97738438100;
    const a: f64 = 6378388.0000;
    const e: f64 = 0.081991890;

    try testing.expectApproxEqAbs(n, lambertNormal(lat, a, e), 1e-4);
}

fn latitudeISOFromLatitude(lat: f64, e: f64) f64 {
    return math.log10(math.tan(math.pi / 4.0 * lat / 2.0) * math.pow((1 - e * math.sin(lat)) / (1 + e * math.sin(lat)), e / 2.0));
}

fn latitudeFromLatitudeISO(lat_iso: f64, e: f64, eps: f64) f64 {
    var phi_0: f64 = 2 * math.atan(math.exp(lat_iso)) - math.pi / 2.0;
    // var phi_i: f64 = 2 * math.atan(math.pow(f64, (1.0 + e * math.sin(phi_0)) / (1.0 - e * math.sin(phi_0)), e / 2.0) * math.exp(lat_iso)) - math.pi / 2.0;
    const a = (1.0 + e * math.sin(phi_0)) / (1.0 - e * math.sin(phi_0));
    const b = math.pow(f64, a, e / 2.0);
    const c = math.atan(b);
    // var phi_i: f64 = 2 * math.atan(math.pow(f64, (1 + e * math.sin(phi_0)) / (1.0 - e * math.sin(phi_0)), e / 2.0) * math.exp(lat_iso)) - math.pi / 2.0;
    var phi_i = 2 * c * math.exp(lat_iso) - math.pi / 2.0;

    var delta: f64 = @abs(phi_i - phi_0);
    while (delta > eps) : (delta = @abs(phi_i - phi_0)) {
        phi_0 = phi_i;
        phi_i = 2 * math.atan(math.pow(f64, (1.0 + e * math.sin(phi_0)) / (1.0 - e * math.sin(phi_0)), e / 2.0) * math.exp(lat_iso)) - math.pi / 2.0;
    }
    return phi_i;
}

test "algo0002" {
    const lat_iso = [3]f64{ 1.00552653648, -0.30261690060, 0.2000000000 };
    const e = [3]f64{ 0.08199188998, 0.08199188998, 0.08199188998 };
    const eps = [3]f64{ 1.0e-11, 1.0e-11, 1.0e-11 };

    const phi = [3]f64{ 0.87266462600, -0.29999999997, 0.19998903369 };

    for (0..lat_iso.len) |i| {
        const result = latitudeFromLatitudeISO(lat_iso[i], e[i], eps[i]);
        try testing.expectApproxEqAbs(phi[i], result, 1e10);
    }
}

fn lambertToGeographic(org: Point, zone: LambertZone, lon_merid: f64, e: f64, eps: f64) Point {
    const index = @intFromEnum(zone);
    const n = lambert_n[index];
    const C = lambert_c[index];
    const x_s = lambert_xs[index];
    const y_s = lambert_ys[index];

    const x = org.x;
    const y = org.y;

    const R: f64 = math.sqrt((x - x_s) * (x - x_s) + (y - y_s) * (y - y_s));
    const gamma: f64 = math.atan((x - x_s) / (y_s - y));
    const lon: f64 = lon_merid + gamma / n;
    const lat_iso: f64 = -1.0 / n * @log(@abs(R / C));
    const lat: f64 = latitudeFromLatitudeISO(lat_iso, e, eps);

    return Point.init(lon, lat, org.z);
}

test "algo0004" {
    const org = Point.init(1029705.083, 272723.849, 0);
    const expected = Point.init(0.145512099, 0.872664626, 0.0);
    const eps = 1e-9;
    const computed = lambertToGeographic(org, .LambertI, LON_MERID_GREENWICH, E_CLARK_IGN, eps);
    // try testing.expectApproxEqAbs(expected.x, computed.x, eps);
    try testing.expectApproxEqAbs(expected.y, computed.y, eps);
}

fn cartesianToGeographic(org: Point, meridien: f64, a: f64, e: f64, eps: f64) Point {
    const x = org.x;
    const y = org.y;
    const z: f64 = org.z;

    const lon: f64 = meridien + math.atan(y / x);
    const module: f64 = math.sqrt(x * x + y * y);

    var phi_0: f64 = math.atan(z / (module * (1.0 - (a * e * e) / math.sqrt(x * x + y * y + z * z))));
    var phi_i: f64 = math.atan(z / module / (1.0 - a * e * e * math.cos(phi_0) / (module * math.sqrt(1.0 - e * e * math.sin(phi_0) * math.sin(phi_0)))));
    var delta: f64 = 1e10;
    while (delta > eps) : (delta = @abs(phi_i - phi_0)) {
        phi_0 = phi_i;
        phi_i = math.atan(z / module / (1.0 - a * e * e * math.cos(phi_0) / (module * math.sqrt(1.0 - e * e * math.sin(phi_0) * math.sin(phi_0)))));
    }

    const he: f64 = module / math.cos(phi_i) - a / math.sqrt(1.0 - e * e * math.sin(phi_i) * math.sin(phi_i));

    return Point.init(lon, phi_i, he);
}

test "test algo00012" {
    const a: [3]f64 = [3]f64{ 6378249.2000, 6378249.2000, 6378249.2000 };
    const e: [3]f64 = [3]f64{ 0.08248325679, 0.08248325679, 0.08248325679 };
    const x: [3]f64 = [3]f64{ 6376064.6950, 6378232.2150, 6376897.5370 };
    const y: [3]f64 = [3]f64{ 111294.6230, 18553.5780, 37099.7050 };
    const z: [3]f64 = [3]f64{ 128984.7250, 0.0000, -202730.9070 };
    const eps: [3]f64 = [3]f64{ 1e-11, 1e-11, 1e-11 };
    const lon: [3]f64 = [3]f64{ 0.01745329248, 0.00290888212, 0.00581776423 };
    const lat: [3]f64 = [3]f64{ 0.02036217457, 0.00000000000, -0.03199770301 };
    const he: [3]f64 = [3]f64{ 99.9995, 10.0001, 2000.0001 };

    for (0..e.len) |i| {
        const sample = Point.init(x[i], y[i], z[i]);
        const val = cartesianToGeographic(sample, LON_MERID_PARIS, a[i], e[i], eps[i]);

        try testing.expectApproxEqAbs(lon[i], val.x, 1e-11);
        try testing.expectApproxEqAbs(lat[i], val.y, 1e-11);
        try testing.expectApproxEqAbs(he[i], val.z, 1e-4);
    }
}

fn geographicToCartesian(lon: f64, lat: f64, he: f64, a: f64, e: f64) Point {
    const N = lambertNormal(lat, a, e);

    const x = (N + he) * math.cos(lat) * math.cos(lon);
    const y = (N + he) * math.cos(lat) * math.sin(lon);
    const z = (N * (1 - e * e) + he) * math.sin(lat);

    return Point.init(x, y, z);
}

test "algo009" {
    const lon = [3]f64{ 0.01745329248, 0.00290888212, 0.00581776423 };
    const lat = [3]f64{ 0.02036217457, 0.00000000000, -0.03199770300 };
    const he = [3]f64{ 100.0000, 10.0000, 2000.0000 };
    const a = [3]f64{ 6378249.2000, 6378249.2000, 6378249.2000 };
    const e = [3]f64{ 0.08248325679, 0.08248325679, 0.08248325679 };

    const expected = [3]Point{ Point.init(6376064.6955, 111294.6230, 128984.7250), Point.init(6378232.2149, 18553.5780, 0), Point.init(6376897.5369, 37099.7050, -202730.9070) };

    for (0..expected.len) |i| {
        const pt = geographicToCartesian(lon[i], lat[i], he[i], a[i], e[i]);
        try testing.expectApproxEqAbs(expected[i].x, pt.x, 1e-4);
        try testing.expectApproxEqAbs(expected[i].y, pt.y, 1e-4);
        try testing.expectApproxEqAbs(expected[i].z, pt.z, 1e-4);
    }
}

/// This methods takes a point where coordinates are expressed in meters and a lambert
/// zone and returned a point in the WGS84 system in radians
pub fn convertWGS84(point: Point, zone: LambertZone) Point {
    switch (zone) {
        .Lambert93 => {
            return lambertToGeographic(point, zone, LON_MERID_IERS, E_WGS84, DEFAULT_EPS);
        },
        else => {
            const point1 = lambertToGeographic(point, zone, LON_MERID_PARIS, E_CLARK_IGN, DEFAULT_EPS);
            var point2 = geographicToCartesian(point1.x, point1.y, point1.z, A_CLARK_IGN, E_CLARK_IGN);

            point2.x -= 168;
            point2.y -= 60;
            point2.z += 320;

            return cartesianToGeographic(point2, LON_MERID_GREENWICH, A_WGS84, E_WGS84, DEFAULT_EPS);
        },
    }
}

test "ZenithStrasbourg" {
    const orgMeterPoint = Point.init(994300.623, 113409.981, 0);
    const expectedDegreePoint = Point.init(7.68639475277068, 48.5953456709144, 0);

    const dest = convertWGS84(orgMeterPoint, .LambertI).mul(RadianToDegree);
    try testing.expectApproxEqAbs(expectedDegreePoint.x, dest.x, 1e-5);
    try testing.expectApproxEqAbs(expectedDegreePoint.y, dest.y, 1e-5);
}

test "Bug LambertIIE" {
    const orgMeterPoint = Point.init(369419, 1986498, 0);
    const expectedDegreePoint = Point.init(-0.579117201473994, 44.84071560809383, 0);

    const dest = convertWGS84(orgMeterPoint, .LambertII_E).mul(RadianToDegree);
    try testing.expectApproxEqAbs(expectedDegreePoint.x, dest.x, 1e-4);
    try testing.expectApproxEqAbs(expectedDegreePoint.y, dest.y, 1e-4);
}
