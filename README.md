# Lambert for Zig

A simple zig module to convert Lambert coordinates to GPS WGS84 coordinates based on the [IGN algorithms and methods](http://geodesie.ign.fr/contenu/fichiers/documentation/algorithmes/notice/NTG_71.pdf)

# Usage

```zig
const lambert = @import("lambert");

const point = lambert.Point.init(994300.623, 113409.981, 0);
const converted = lambert.convertWGS84(orgMeterPoint, .LambertI).mul(lambert.RadianToDegree);

```