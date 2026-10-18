/*
* Copyright (c) 2006-2007 Erin Catto http://www.gphysics.com
*
* This software is provided 'as-is', without any express or implied
* warranty.  In no event will the authors be held liable for any damages
* arising from the use of this software.
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely, subject to the following restrictions:
* 1. The origin of this software must not be misrepresented; you must not
* claim that you wrote the original software. If you use this software
* in a product, an acknowledgment in the product documentation would be
* appreciated but is not required.
* 2. Altered source versions must be plainly marked as such, and must not be
* misrepresented as being the original software.
* 3. This notice may not be removed or altered from any source distribution.
*/

const b2Color = class {
  constructor(rr = 0, gg = 0, bb = 0) {
    this._r = Math.trunc(255 * b2Math.Clamp(rr, 0.0, 1.0));
    this._g = Math.trunc(255 * b2Math.Clamp(gg, 0.0, 1.0));
    this._b = Math.trunc(255 * b2Math.Clamp(bb, 0.0, 1.0));
  }
  Set(rr = 0, gg = 0, bb = 0) {
    this._r = Math.trunc(255 * b2Math.Clamp(rr, 0.0, 1.0));
    this._g = Math.trunc(255 * b2Math.Clamp(gg, 0.0, 1.0));
    this._b = Math.trunc(255 * b2Math.Clamp(bb, 0.0, 1.0));
  }
  set r(rr = 0) {
    this._r = Math.trunc(255 * b2Math.Clamp(rr, 0.0, 1.0));
  }
  set g(gg = 0) {
    this._g = Math.trunc(255 * b2Math.Clamp(gg, 0.0, 1.0));
  }
  set b(bb = 0) {
    this._b = Math.trunc(255 * b2Math.Clamp(bb, 0.0, 1.0));
  }
  get color() {
    return (this._r << 16) | (this._g << 8) | this._b;
  }
};

const b2Settings = class {
  static VERSION = "2.1alpha";
  static USHRT_MAX = 0x0000ffff;
  static b2_pi = Math.PI;
  static b2_maxManifoldPoints = 2;
  static b2_aabbExtension = 0.1;
  static b2_aabbMultiplier = 2.0;
  static b2_epsilon = 1.1920928955078125e-7;
  static b2_linearSlop = 0.005;
  static b2_polygonRadius = 2.0 * this.b2_linearSlop;
  static b2_angularSlop = 2.0 / 180.0 * this.b2_pi;
  static b2_toiSlop = 8.0 * this.b2_linearSlop;
  static b2_maxTOIContactsPerIsland = 32;
  static b2_maxTOIJointsPerIsland = 32;
  static b2_velocityThreshold = 1.0;
  static b2_maxLinearCorrection = 0.2;
  static b2_maxAngularCorrection = 8.0 / 180.0 * this.b2_pi;
  static b2_maxTranslation = 2.0;
  static b2_maxTranslationSquared = this.b2_maxTranslation * this.b2_maxTranslation;
  static b2_maxRotation = 0.5 * this.b2_pi;
  static b2_maxRotationSquared = this.b2_maxRotation * this.b2_maxRotation;
  static b2_contactBaumgarte = 0.2;
  static b2_timeToSleep = 0.5;
  static b2_linearSleepTolerance = 0.01;
  static b2_angularSleepTolerance = 2.0 / 180.0 * this.b2_pi;
  static b2MixFriction(friction1 = 0, friction2 = 0) {
    return Math.sqrt(friction1 * friction2);
  }
  static b2MixRestitution(restitution1 = 0, restitution2 = 0) {
    return restitution1 > restitution2 ? restitution1 : restitution2;
  }
  static b2Assert(a) {
    if (!a) {
      throw new Error("Assertion Failed");
    }
  }
};

const b2Mat22 = class {
  constructor() {
    this.col1 = new b2Vec2();
    this.col2 = new b2Vec2();
    this.SetIdentity();
  }
  static FromAngle(angle = 0) {
    const mat = new b2Mat22();
    mat.Set(angle);
    return mat;
  }
  static FromVV(c1, c2) {
    const mat = new b2Mat22();
    mat.SetVV(c1, c2);
    return mat;
  }
  Set(angle = 0) {
    const sinCos = b2Math.s_sinCos;
    b2Math.SinCos(angle, sinCos);
    const c = sinCos[1];
    const s = sinCos[0];
    this.col1.x = c;
    this.col2.x = -s;
    this.col1.y = s;
    this.col2.y = c;
  }
  SetVV(c1, c2) {
    this.col1.SetV(c1);
    this.col2.SetV(c2);
  }
  Copy() {
    const mat = new b2Mat22();
    mat.SetM(this);
    return mat;
  }
  SetM(m) {
    this.col1.SetV(m.col1);
    this.col2.SetV(m.col2);
  }
  AddM(m) {
    this.col1.x += m.col1.x;
    this.col1.y += m.col1.y;
    this.col2.x += m.col2.x;
    this.col2.y += m.col2.y;
  }
  SetIdentity() {
    this.col1.x = 1.0;
    this.col2.x = 0.0;
    this.col1.y = 0.0;
    this.col2.y = 1.0;
  }
  SetZero() {
    this.col1.x = 0.0;
    this.col2.x = 0.0;
    this.col1.y = 0.0;
    this.col2.y = 0.0;
  }
  GetAngle() {
    return b2Math.Atan2(this.col1.y, this.col1.x);
  }
  GetInverse(out) {
    const a = this.col1.x;
    const b = this.col2.x;
    const c = this.col1.y;
    const d = this.col2.y;
    let det = a * d - b * c;
    if (det != 0.0) {
      det = 1.0 / det;
    }
    out.col1.x = det * d;
    out.col2.x = -det * b;
    out.col1.y = -det * c;
    out.col2.y = det * a;
    return out;
  }
  Solve(out, bX = 0, bY = 0) {
    const a11 = this.col1.x;
    const a12 = this.col2.x;
    const a21 = this.col1.y;
    const a22 = this.col2.y;
    let det = a11 * a22 - a12 * a21;
    if (det != 0.0) {
      det = 1.0 / det;
    }
    out.x = det * (a22 * bX - a12 * bY);
    out.y = det * (a11 * bY - a21 * bX);
    return out;
  }
  Abs() {
    this.col1.Abs();
    this.col2.Abs();
  }
};

const b2Mat33 = class {
  constructor(c1 = null, c2 = null, c3 = null) {
    this.col1 = new b2Vec3();
    this.col2 = new b2Vec3();
    this.col3 = new b2Vec3();
    if (!c1 && !c2 && !c3) {
      this.col1.SetZero();
      this.col2.SetZero();
      this.col3.SetZero();
    } else {
      this.col1.SetV(c1);
      this.col2.SetV(c2);
      this.col3.SetV(c3);
    }
  }
  SetVVV(c1, c2, c3) {
    this.col1.SetV(c1);
    this.col2.SetV(c2);
    this.col3.SetV(c3);
  }
  Copy() {
    return new b2Mat33(this.col1, this.col2, this.col3);
  }
  SetM(m) {
    this.col1.SetV(m.col1);
    this.col2.SetV(m.col2);
    this.col3.SetV(m.col3);
  }
  AddM(m) {
    this.col1.x += m.col1.x;
    this.col1.y += m.col1.y;
    this.col1.z += m.col1.z;
    this.col2.x += m.col2.x;
    this.col2.y += m.col2.y;
    this.col2.z += m.col2.z;
    this.col3.x += m.col3.x;
    this.col3.y += m.col3.y;
    this.col3.z += m.col3.z;
  }
  SetIdentity() {
    this.col1.x = 1.0;
    this.col2.x = 0.0;
    this.col3.x = 0.0;
    this.col1.y = 0.0;
    this.col2.y = 1.0;
    this.col3.y = 0.0;
    this.col1.z = 0.0;
    this.col2.z = 0.0;
    this.col3.z = 1.0;
  }
  SetZero() {
    this.col1.x = 0.0;
    this.col2.x = 0.0;
    this.col3.x = 0.0;
    this.col1.y = 0.0;
    this.col2.y = 0.0;
    this.col3.y = 0.0;
    this.col1.z = 0.0;
    this.col2.z = 0.0;
    this.col3.z = 0.0;
  }
  Solve22(out, bX = 0, bY = 0) {
    const a11 = this.col1.x;
    const a12 = this.col2.x;
    const a21 = this.col1.y;
    const a22 = this.col2.y;
    let det = a11 * a22 - a12 * a21;
    if (det != 0.0) {
      det = 1.0 / det;
    }
    out.x = det * (a22 * bX - a12 * bY);
    out.y = det * (a11 * bY - a21 * bX);
    return out;
  }
  Solve33(out, bX = 0, bY = 0, bZ = 0) {
    const a11 = this.col1.x;
    const a21 = this.col1.y;
    const a31 = this.col1.z;
    const a12 = this.col2.x;
    const a22 = this.col2.y;
    const a32 = this.col2.z;
    const a13 = this.col3.x;
    const a23 = this.col3.y;
    const a33 = this.col3.z;
    let det = a11 * (a22 * a33 - a32 * a23) + a21 * (a32 * a13 - a12 * a33) + a31 * (a12 * a23 - a22 * a13);
    if (det != 0.0) {
      det = 1.0 / det;
    }
    out.x = det * (bX * (a22 * a33 - a32 * a23) + bY * (a32 * a13 - a12 * a33) + bZ * (a12 * a23 - a22 * a13));
    out.y = det * (a11 * (bY * a33 - bZ * a23) + a21 * (bZ * a13 - bX * a33) + a31 * (bX * a23 - bY * a13));
    out.z = det * (a11 * (a22 * bZ - a32 * bY) + a21 * (a32 * bX - a12 * bZ) + a31 * (a12 * bY - a22 * bX));
    return out;
  }
};

const b2Vec2 = class {
  constructor(x_ = 0, y_ = 0) {
    this.x = x_;
    this.y = y_;
  }
  SetZero() {
    this.x = 0.0;
    this.y = 0.0;
  }
  Set(x_ = 0, y_ = 0) {
    this.x = x_;
    this.y = y_;
  }
  SetV(v) {
    this.x = v.x;
    this.y = v.y;
  }
  GetNegative() {
    return new b2Vec2(-this.x, -this.y);
  }
  NegativeSelf() {
    this.x = -this.x;
    this.y = -this.y;
  }
  static Make(x_ = 0, y_ = 0) {
    return new b2Vec2(x_, y_);
  }
  Copy() {
    return new b2Vec2(this.x, this.y);
  }
  Add(v) {
    this.x += v.x;
    this.y += v.y;
  }
  Subtract(v) {
    this.x -= v.x;
    this.y -= v.y;
  }
  Multiply(a = 0) {
    this.x *= a;
    this.y *= a;
  }
  MulM(A) {
    const tX = this.x;
    this.x = A.col1.x * tX + A.col2.x * this.y;
    this.y = A.col1.y * tX + A.col2.y * this.y;
  }
  MulTM(A) {
    const tX = b2Math.Dot(this, A.col1);
    this.y = b2Math.Dot(this, A.col2);
    this.x = tX;
  }
  CrossVF(s = 0) {
    const tX = this.x;
    this.x = s * this.y;
    this.y = -s * tX;
  }
  CrossFV(s = 0) {
    const tX = this.x;
    this.x = -s * this.y;
    this.y = s * tX;
  }
  MinV(b) {
    this.x = this.x < b.x ? this.x : b.x;
    this.y = this.y < b.y ? this.y : b.y;
  }
  MaxV(b) {
    this.x = this.x > b.x ? this.x : b.x;
    this.y = this.y > b.y ? this.y : b.y;
  }
  Abs() {
    if (this.x < 0) this.x = -this.x;
    if (this.y < 0) this.y = -this.y;
  }
  Length() {
    return Math.sqrt(this.x * this.x + this.y * this.y);
  }
  LengthSquared() {
    return this.x * this.x + this.y * this.y;
  }
  Normalize() {
    const length = Math.sqrt(this.x * this.x + this.y * this.y);
    if (length < b2Settings.b2_epsilon) {
      return 0.0;
    }
    const invLength = 1.0 / length;
    this.x *= invLength;
    this.y *= invLength;
    return length;
  }
  IsValid() {
    return b2Math.IsValid(this.x) && b2Math.IsValid(this.y);
  }
};

const b2Transform = class {
  constructor(pos = null, r = null) {
    this.position = new b2Vec2();
    this.R = new b2Mat22();
    if (pos) {
      this.position.SetV(pos);
      this.R.SetM(r);
    }
  }
  Initialize(pos, r) {
    this.position.SetV(pos);
    this.R.SetM(r);
  }
  SetIdentity() {
    this.position.SetZero();
    this.R.SetIdentity();
  }
  Set(x) {
    this.position.SetV(x.position);
    this.R.SetM(x.R);
  }
  GetAngle() {
    return b2Math.Atan2(this.R.col1.y, this.R.col1.x);
  }
};

const b2Math = class {
  static s_reduction = new Float64Array(2);
  static twoOverPiBits = new Float64Array([
    0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
    0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
    0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
    0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
    0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
    0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
    0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
    0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
    0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
    0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
    0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B
  ]);
  static pio2Parts = new Float64Array([1.570796251296997, 7.549789415861596e-8, 5.390302529957765e-15, 3.282003415807913e-22, 1.270655753080676e-29, 1.2293330898111133e-36, 2.7337005381646456e-44, 2.1674168387780482e-51]);
  static atanHigh = new Float64Array([0.4636476090008061, 0.7853981633974483, 0.982793723247329, 1.5707963267948966]);
  static atanLow = new Float64Array([2.2698777452961687e-17, 3.061616997868383e-17, 1.3903311031230998e-17, 6.123233995736766e-17]);
  static s_angleChunks = new Float64Array(3);
  static s_reductionTerms = new Float64Array(20);
  static s_reductionProducts = new Float64Array(20);
  static s_reductionPiProducts = new Float64Array(20);
  static s_reductionProductChunks = new Int32Array(20);
  static s_sinCos = [0.0, 0.0];
  static s_float64 = new Float64Array(1);
  static s_words = new Int32Array(this.s_float64.buffer);
  static s_highWordIndex = new Uint8Array(new Uint16Array([1]).buffer)[0] == 1 ? 1 : 0;
  static powersOfTwo = new Float64Array(2046);
  static {
    this.powersOfTwo[0] = 2.2250738585072014e-308;
    for (let i = 1; i < 2046; ++i) {
      this.powersOfTwo[i] = this.powersOfTwo[i - 1] * 2.0;
    }
  }
  static b2Vec2_zero = new b2Vec2(0.0, 0.0);
  static b2Mat22_identity = b2Mat22.FromVV(new b2Vec2(1.0, 0.0), new b2Vec2(0.0, 1.0));
  static b2Transform_identity = new b2Transform(this.b2Vec2_zero, this.b2Mat22_identity);
  static IsValid(x = 0) {
    return isFinite(x);
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y;
  }
  static CrossVV(a, b) {
    return a.x * b.y - a.y * b.x;
  }
  static CrossVF(a, s = 0) {
    const v = new b2Vec2(s * a.y, (-s * a.x));
    return v;
  }
  static CrossFV(s = 0, a) {
    const v = new b2Vec2((-s * a.y), s * a.x);
    return v;
  }
  static MulMV(A, v) {
    const u = new b2Vec2(A.col1.x * v.x + A.col2.x * v.y, A.col1.y * v.x + A.col2.y * v.y);
    return u;
  }
  static MulTMV(A, v) {
    const u = new b2Vec2(b2Math.Dot(v, A.col1), b2Math.Dot(v, A.col2));
    return u;
  }
  static MulX(T, v) {
    const a = b2Math.MulMV(T.R, v);
    a.x += T.position.x;
    a.y += T.position.y;
    return a;
  }
  static MulXT(T, v) {
    const a = b2Math.SubtractVV(v, T.position);
    const tX = a.x * T.R.col1.x + a.y * T.R.col1.y;
    a.y = a.x * T.R.col2.x + a.y * T.R.col2.y;
    a.x = tX;
    return a;
  }
  static AddVV(a, b) {
    const v = new b2Vec2(a.x + b.x, a.y + b.y);
    return v;
  }
  static SubtractVV(a, b) {
    const v = new b2Vec2(a.x - b.x, a.y - b.y);
    return v;
  }
  static Distance(a, b) {
    const cX = a.x - b.x;
    const cY = a.y - b.y;
    return Math.sqrt(cX * cX + cY * cY);
  }
  static DistanceSquared(a, b) {
    const cX = a.x - b.x;
    const cY = a.y - b.y;
    return cX * cX + cY * cY;
  }
  static MulFV(s = 0, a) {
    const v = new b2Vec2(s * a.x, s * a.y);
    return v;
  }
  static AddMM(A, B) {
    const C = b2Mat22.FromVV(b2Math.AddVV(A.col1, B.col1), b2Math.AddVV(A.col2, B.col2));
    return C;
  }
  static MulMM(A, B) {
    const C = b2Mat22.FromVV(b2Math.MulMV(A, B.col1), b2Math.MulMV(A, B.col2));
    return C;
  }
  static MulTMM(A, B) {
    const c1 = new b2Vec2(b2Math.Dot(A.col1, B.col1), b2Math.Dot(A.col2, B.col1));
    const c2 = new b2Vec2(b2Math.Dot(A.col1, B.col2), b2Math.Dot(A.col2, B.col2));
    const C = b2Mat22.FromVV(c1, c2);
    return C;
  }
  static Abs(a = 0) {
    return a > 0.0 ? a : -a;
  }
  static AbsV(a) {
    const b = new b2Vec2(b2Math.Abs(a.x), b2Math.Abs(a.y));
    return b;
  }
  static AbsM(A) {
    const B = b2Mat22.FromVV(b2Math.AbsV(A.col1), b2Math.AbsV(A.col2));
    return B;
  }
  static Min(a = 0, b = 0) {
    return a < b ? a : b;
  }
  static MinV(a, b) {
    const c = new b2Vec2(b2Math.Min(a.x, b.x), b2Math.Min(a.y, b.y));
    return c;
  }
  static Max(a = 0, b = 0) {
    return a > b ? a : b;
  }
  static MaxV(a, b) {
    const c = new b2Vec2(b2Math.Max(a.x, b.x), b2Math.Max(a.y, b.y));
    return c;
  }
  static Clamp(a = 0, low = 0, high = 0) {
    return a < low ? low : a > high ? high : a;
  }
  static ClampV(a, low, high) {
    return b2Math.MaxV(low, b2Math.MinV(a, high));
  }
  static Swap(a, b) {
    const tmp = a[0];
    a[0] = b[0];
    b[0] = tmp;
  }
  static Random() {
    return Math.random() * 2 - 1;
  }
  static RandomRange(lo = 0, hi = 0) {
    let r = Math.random();
    r = (hi - lo) * r + lo;
    return r;
  }
  static NextPowerOfTwo(x = 0) {
    x |= (x >> 1) & 0x7FFFFFFF;
    x |= (x >> 2) & 0x3FFFFFFF;
    x |= (x >> 4) & 0x0FFFFFFF;
    x |= (x >> 8) & 0x00FFFFFF;
    x |= (x >> 16) & 0x0000FFFF;
    return x + 1;
  }
  static IsPowerOfTwo(x = 0) {
    const result = x > 0 && (x & (x - 1)) == 0;
    return result;
  }
  /*
   * JavaScript ports of fdlibm's sin, cos and atan2 to secure browser engine parity.
   *
   * ====================================================
   * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
   *
   * Developed at SunSoft, a Sun Microsystems, Inc. business.
   * Permission to use, copy, modify, and distribute this
   * software is freely granted, provided that this notice
   * is preserved.
   * ====================================================
   */
  static Sin(x = 0) {
    if (Math.abs(x) < 1.4901161193847656e-8) return x;
    if (Math.abs(x) < 0.7853984832763672) return b2Math.SinKernel(x, 0.0, false);
    if (!b2Math.IsValid(x)) return x - x;
    const reduction = b2Math.s_reduction;
    switch (b2Math.ReduceAngle(x, reduction) & 3) {
    case 0:
      return b2Math.SinKernel(reduction[0], reduction[1], true);
    case 1:
      return b2Math.CosKernel(reduction[0], reduction[1]);
    case 2:
      return -b2Math.SinKernel(reduction[0], reduction[1], true);
    default:
      return -b2Math.CosKernel(reduction[0], reduction[1]);
    }
  }
  static Cos(x = 0) {
    if (Math.abs(x) < 0.7853984832763672) return b2Math.CosKernel(x, 0.0);
    if (!b2Math.IsValid(x)) return x - x;
    const reduction = b2Math.s_reduction;
    switch (b2Math.ReduceAngle(x, reduction) & 3) {
    case 0:
      return b2Math.CosKernel(reduction[0], reduction[1]);
    case 1:
      return -b2Math.SinKernel(reduction[0], reduction[1], true);
    case 2:
      return -b2Math.CosKernel(reduction[0], reduction[1]);
    default:
      return b2Math.SinKernel(reduction[0], reduction[1], true);
    }
  }
  static SinCos(x = 0, sinCos) {
    if (Math.abs(x) < 0.7853984832763672) {
      sinCos[0] = Math.abs(x) < 1.4901161193847656e-8 ? x : b2Math.SinKernel(x, 0.0, false);
      sinCos[1] = b2Math.CosKernel(x, 0.0);
      return;
    }
    if (!b2Math.IsValid(x)) {
      sinCos[0] = x - x;
      sinCos[1] = x - x;
      return;
    }
    const reduction = b2Math.s_reduction;
    const n = b2Math.ReduceAngle(x, reduction);
    const s = b2Math.SinKernel(reduction[0], reduction[1], true);
    const c = b2Math.CosKernel(reduction[0], reduction[1]);
    switch (n & 3) {
    case 0:
      sinCos[0] = s;
      sinCos[1] = c;
      break;
    case 1:
      sinCos[0] = c;
      sinCos[1] = -s;
      break;
    case 2:
      sinCos[0] = -s;
      sinCos[1] = -c;
      break;
    default:
      sinCos[0] = -c;
      sinCos[1] = s;
      break;
    }
  }
  static ReduceAngle(x = 0, reduction) {
    const k_invPio2 = 0.6366197723675814;
    const k_pio2Part1 = 1.5707963267341256;
    const k_pio2Tail1 = 6.077100506506192e-11;
    const k_pio2Part2 = 6.077100506303966e-11;
    const k_pio2Tail2 = 2.0222662487959506e-21;
    const k_pio2Part3 = 2.0222662487111665e-21;
    const k_pio2Tail3 = 8.4784276603689e-32;
    const k_roundingShift = 6755399441055744.0;
    const high = b2Math.HighWord(x);
    const absHigh = high & 0x7FFFFFFF;
    let n = 0;
    let z;
    if (absHigh <= 0x401C463B) {
      if (absHigh <= 0x400F6A7A) {
        if ((absHigh & 0xFFFFF) != 0x921FB) n = absHigh <= 0x4002D97C ? 1 : 2;
      } else if (absHigh <= 0x4015FDBC) {
        if (absHigh != 0x4012D97C) n = 3;
      } else if (absHigh != 0x401921FB) {
        n = 4;
      }
      if (n != 0) {
        if (high > 0) {
          z = x - n * k_pio2Part1;
          reduction[0] = z - n * k_pio2Tail1;
          reduction[1] = (z - reduction[0]) - n * k_pio2Tail1;
          return n;
        }
        z = x + n * k_pio2Part1;
        reduction[0] = z + n * k_pio2Tail1;
        reduction[1] = (z - reduction[0]) + n * k_pio2Tail1;
        return -n;
      }
    }
    if (absHigh >= 0x413921FB) return b2Math.ReduceLargeAngle(x, reduction);
    n = (x * k_invPio2 + k_roundingShift) - k_roundingShift;
    const exponent = absHigh >> 20;
    let r = x - n * k_pio2Part1;
    let w = n * k_pio2Tail1;
    let t;
    reduction[0] = r - w;
    if (exponent - ((b2Math.HighWord(reduction[0]) >> 20) & 0x7FF) > 16) {
      t = r;
      w = n * k_pio2Part2;
      r = t - w;
      w = n * k_pio2Tail2 - ((t - r) - w);
      reduction[0] = r - w;
      if (exponent - ((b2Math.HighWord(reduction[0]) >> 20) & 0x7FF) > 49) {
        t = r;
        w = n * k_pio2Part3;
        r = t - w;
        w = n * k_pio2Tail3 - ((t - r) - w);
        reduction[0] = r - w;
      }
    }
    reduction[1] = (r - reduction[0]) - w;
    return n;
  }
  static ReduceLargeAngle(x = 0, reduction) {
    const k_two24 = 16777216.0;
    const chunks = b2Math.s_angleChunks;
    const exponent = ((b2Math.HighWord(x) & 0x7FFFFFFF) >> 20) - 1046;
    let z = Math.abs(x) * b2Math.PowerOfTwo(-exponent);
    chunks[0] = Math.floor(z);
    z = (z - chunks[0]) * k_two24;
    chunks[1] = Math.floor(z);
    chunks[2] = (z - chunks[1]) * k_two24;
    let chunkCount = 3;
    while (chunks[chunkCount - 1] == 0.0) {
      --chunkCount;
    }
    const n = b2Math.ReduceLargeAngleKernel(chunks, chunkCount, exponent, reduction);
    if (x < 0.0) {
      reduction[0] = -reduction[0];
      reduction[1] = -reduction[1];
      return -n;
    }
    return n;
  }
  static ReduceLargeAngleKernel(chunks, chunkCount = 0, exponent = 0, reduction) {
    const k_two24 = 16777216.0;
    const k_twoNegative24 = 5.960464477539063e-8;
    const k_initialTerms = 4;
    let i;
    let j;
    let k;
    let n;
    let half;
    let carry;
    let sum;
    let high;
    let z;
    const k_twoOverPi = b2Math.twoOverPiBits;
    const k_pio2 = b2Math.pio2Parts;
    const terms = b2Math.s_reductionTerms;
    const products = b2Math.s_reductionProducts;
    const piProducts = b2Math.s_reductionPiProducts;
    const productChunks = b2Math.s_reductionProductChunks;
    const lastChunk = chunkCount - 1;
    let lastProduct = k_initialTerms;
    let tableStart = Math.floor((exponent - 3) / 24);
    if (tableStart < 0) tableStart = 0;
    let scale = exponent - 24 * (tableStart + 1);
    for (i = 0; i <= lastChunk + k_initialTerms; ++i) {
      j = tableStart - lastChunk + i;
      terms[i] = j < 0 ? 0.0 : k_twoOverPi[j];
    }
    for (i = 0; i <= k_initialTerms; ++i) {
      sum = 0.0;
      for (j = 0; j <= lastChunk; ++j) {
        sum += chunks[j] * terms[lastChunk + i - j];
      }
      products[i] = sum;
    }
    for (;;) {
      z = products[lastProduct];
      for (i = 0; i < lastProduct; ++i) {
        high = Math.floor(k_twoNegative24 * z);
        productChunks[i] = Math.floor(z - k_two24 * high);
        z = products[lastProduct - i - 1] + high;
      }
      z *= b2Math.PowerOfTwo(scale);
      z -= 8.0 * Math.floor(z * 0.125);
      n = Math.floor(z);
      z -= n;
      half = 0;
      if (scale > 0) {
        i = productChunks[lastProduct - 1] >> (24 - scale);
        n += i;
        productChunks[lastProduct - 1] -= i << (24 - scale);
        half = productChunks[lastProduct - 1] >> (23 - scale);
      } else if (scale == 0) {
        half = productChunks[lastProduct - 1] >> 23;
      } else if (z >= 0.5) {
        half = 2;
      }
      if (half > 0) {
        n += 1;
        carry = 0;
        for (i = 0; i < lastProduct; ++i) {
          j = productChunks[i];
          if (carry == 0) {
            if (j != 0) {
              carry = 1;
              productChunks[i] = 0x1000000 - j;
            }
          } else {
            productChunks[i] = 0xFFFFFF - j;
          }
        }
        switch (scale) {
        case 1:
          productChunks[lastProduct - 1] &= 0x7FFFFF;
          break;
        case 2:
          productChunks[lastProduct - 1] &= 0x3FFFFF;
          break;
        }
        if (half == 2) {
          z = 1.0 - z;
          if (carry != 0) z -= b2Math.PowerOfTwo(scale);
        }
      }
      if (z != 0.0) break;
      j = 0;
      for (i = lastProduct - 1; i >= k_initialTerms; --i) {
        j |= productChunks[i];
      }
      if (j != 0) break;
      k = 1;
      while (productChunks[k_initialTerms - k] == 0) {
        ++k;
      }
      for (i = lastProduct + 1; i <= lastProduct + k; ++i) {
        terms[lastChunk + i] = k_twoOverPi[tableStart + i];
        sum = 0.0;
        for (j = 0; j <= lastChunk; ++j) {
          sum += chunks[j] * terms[lastChunk + i - j];
        }
        products[i] = sum;
      }
      lastProduct += k;
    }
    if (z == 0.0) {
      lastProduct -= 1;
      scale -= 24;
      while (productChunks[lastProduct] == 0) {
        lastProduct -= 1;
        scale -= 24;
      }
    } else {
      z *= b2Math.PowerOfTwo(-scale);
      if (z >= k_two24) {
        high = Math.floor(k_twoNegative24 * z);
        productChunks[lastProduct] = Math.floor(z - k_two24 * high);
        lastProduct += 1;
        scale += 24;
        productChunks[lastProduct] = high;
      } else {
        productChunks[lastProduct] = Math.floor(z);
      }
    }
    let power = b2Math.PowerOfTwo(scale);
    for (i = lastProduct; i >= 0; --i) {
      products[i] = power * productChunks[i];
      power *= k_twoNegative24;
    }
    for (i = lastProduct; i >= 0; --i) {
      sum = 0.0;
      for (k = 0; k <= k_initialTerms && k <= lastProduct - i; ++k) {
        sum += k_pio2[k] * products[i + k];
      }
      piProducts[lastProduct - i] = sum;
    }
    sum = 0.0;
    for (i = lastProduct; i >= 0; --i) {
      sum += piProducts[i];
    }
    reduction[0] = half == 0 ? sum : -sum;
    sum = piProducts[0] - sum;
    for (i = 1; i <= lastProduct; ++i) {
      sum += piProducts[i];
    }
    reduction[1] = half == 0 ? sum : -sum;
    return n & 7;
  }
  static SinKernel(x = 0, tail = 0, hasTail) {
    const k_sin1 = -0.16666666666666632;
    const k_sin2 = 0.00833333333332249;
    const k_sin3 = -0.0001984126982985795;
    const k_sin4 = 0.0000027557313707070068;
    const k_sin5 = -2.5050760253406863e-8;
    const k_sin6 = 1.58969099521155e-10;
    const z = x * x;
    const w = z * z;
    const r = k_sin2 + z * (k_sin3 + z * k_sin4) + z * w * (k_sin5 + z * k_sin6);
    const v = z * x;
    if (!hasTail) return x + v * (k_sin1 + z * r);
    return x - ((z * (0.5 * tail - v * r) - tail) - v * k_sin1);
  }
  static CosKernel(x = 0, tail = 0) {
    const k_cos1 = 0.0416666666666666;
    const k_cos2 = -0.001388888888887411;
    const k_cos3 = 0.00002480158728947673;
    const k_cos4 = -2.7557314351390663e-7;
    const k_cos5 = 2.087572321298175e-9;
    const k_cos6 = -1.1359647557788195e-11;
    const z = x * x;
    let w = z * z;
    const r = z * (k_cos1 + z * (k_cos2 + z * k_cos3)) + w * w * (k_cos4 + z * (k_cos5 + z * k_cos6));
    const hz = 0.5 * z;
    w = 1.0 - hz;
    return w + (((1.0 - w) - hz) + (z * r - x * tail));
  }
  static Atan2(y = 0, x = 0) {
    const k_pi = 3.141592653589793;
    const k_piOver2 = 1.5707963267948966;
    const k_piOver4 = 0.7853981633974483;
    const k_piLow = 1.2246467991473532e-16;
    if (isNaN(x) || isNaN(y)) return x + y;
    if (x == 1.0) return b2Math.Atan(y);
    const highX = b2Math.HighWord(x);
    const highY = b2Math.HighWord(y);
    let signs = ((highY >> 31) & 1) | ((highX >> 30) & 2);
    if (y == 0.0) {
      switch (signs) {
      case 0:
      case 1:
        return y;
      case 2:
        return k_pi;
      default:
        return -k_pi;
      }
    }
    if (x == 0.0) return highY < 0 ? -k_piOver2 : k_piOver2;
    if (!b2Math.IsValid(x)) {
      if (!b2Math.IsValid(y)) {
        switch (signs) {
        case 0:
          return k_piOver4;
        case 1:
          return -k_piOver4;
        case 2:
          return 3.0 * k_piOver4;
        default:
          return -3.0 * k_piOver4;
        }
      }
      switch (signs) {
      case 0:
        return 0.0;
      case 1:
        return -0.0;
      case 2:
        return k_pi;
      default:
        return -k_pi;
      }
    }
    if (!b2Math.IsValid(y)) return highY < 0 ? -k_piOver2 : k_piOver2;
    const exponentDifference = ((highY & 0x7FFFFFFF) - (highX & 0x7FFFFFFF)) >> 20;
    let angle;
    if (exponentDifference > 60) {
      angle = k_piOver2 + 0.5 * k_piLow;
      signs &= 1;
    } else if (highX < 0 && exponentDifference < -60) {
      angle = 0.0;
    } else {
      angle = b2Math.Atan(Math.abs(y / x));
    }
    switch (signs) {
    case 0:
      return angle;
    case 1:
      return -angle;
    case 2:
      return k_pi - (angle - k_piLow);
    default:
      return (angle - k_piLow) - k_pi;
    }
  }
  static Atan(x = 0) {
    const k_atanHigh = b2Math.atanHigh;
    const k_atanLow = b2Math.atanLow;
    const k_atan1 = 0.3333333333333293;
    const k_atan2 = -0.19999999999876483;
    const k_atan3 = 0.14285714272503466;
    const k_atan4 = -0.11111110405462356;
    const k_atan5 = 0.09090887133436507;
    const k_atan6 = -0.0769187620504483;
    const k_atan7 = 0.06661073137387531;
    const k_atan8 = -0.058335701337905735;
    const k_atan9 = 0.049768779946159324;
    const k_atan10 = -0.036531572744216916;
    const k_atan11 = 0.016285820115365782;
    if (isNaN(x)) return x + x;
    const negative = x < 0.0;
    let range = -1;
    if (Math.abs(x) >= 73786976294838210000) return negative ? (-k_atanHigh[3] - k_atanLow[3]) : k_atanHigh[3] + k_atanLow[3];
    if (Math.abs(x) < 7.450580596923828e-9) return x;
    if (Math.abs(x) >= 0.4375) {
      x = Math.abs(x);
      if (x < 1.1875) {
        if (x < 0.6875) {
          range = 0;
          x = (2.0 * x - 1.0) / (2.0 + x);
        } else {
          range = 1;
          x = (x - 1.0) / (x + 1.0);
        }
      } else {
        if (x < 2.4375) {
          range = 2;
          x = (x - 1.5) / (1.0 + 1.5 * x);
        } else {
          range = 3;
          x = -1.0 / x;
        }
      }
    }
    let z = x * x;
    const w = z * z;
    const s1 = z * (k_atan1 + w * (k_atan3 + w * (k_atan5 + w * (k_atan7 + w * (k_atan9 + w * k_atan11)))));
    const s2 = w * (k_atan2 + w * (k_atan4 + w * (k_atan6 + w * (k_atan8 + w * k_atan10))));
    if (range < 0) return x - x * (s1 + s2);
    z = k_atanHigh[range] - ((x * (s1 + s2) - k_atanLow[range]) - x);
    return negative ? -z : z;
  }
  static HighWord(x = 0) {
    b2Math.s_float64[0] = x;
    return b2Math.s_words[b2Math.s_highWordIndex];
  }
  static PowerOfTwo(exponent = 0) {
    if (exponent >= -1022 && exponent <= 1023) return b2Math.powersOfTwo[(exponent + 1022) | 0];
    b2Math.s_words[b2Math.s_highWordIndex] = (exponent + 1023) << 20;
    b2Math.s_words[1 - b2Math.s_highWordIndex] = 0;
    return b2Math.s_float64[0];
  }
};

const b2Sweep = class {
  constructor() {
    this.localCenter = new b2Vec2();
    this.c0 = new b2Vec2();
    this.c = new b2Vec2();
    this.a0 = 0.0;
    this.a = 0.0;
    this.t0 = 0.0;
  }
  Set(other) {
    this.localCenter.SetV(other.localCenter);
    this.c0.SetV(other.c0);
    this.c.SetV(other.c);
    this.a0 = other.a0;
    this.a = other.a;
    this.t0 = other.t0;
  }
  Copy() {
    const copy = new b2Sweep();
    copy.localCenter.SetV(this.localCenter);
    copy.c0.SetV(this.c0);
    copy.c.SetV(this.c);
    copy.a0 = this.a0;
    copy.a = this.a;
    copy.t0 = this.t0;
    return copy;
  }
  GetTransform(xf, alpha = 0) {
    const oneMinusAlpha = 1.0 - alpha;
    const tMat = xf.R;
    tMat.Set(oneMinusAlpha * this.a0 + alpha * this.a);
    const tVec = this.localCenter;
    xf.position.x = oneMinusAlpha * this.c0.x + alpha * this.c.x - (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    xf.position.y = oneMinusAlpha * this.c0.y + alpha * this.c.y - (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
  }
  Advance(t = 0) {
    if (this.t0 < t && 1.0 - this.t0 > b2Settings.b2_epsilon) {
      const alpha = (t - this.t0) / (1.0 - this.t0);
      const oneMinusAlpha = 1.0 - alpha;
      this.c0.x = oneMinusAlpha * this.c0.x + alpha * this.c.x;
      this.c0.y = oneMinusAlpha * this.c0.y + alpha * this.c.y;
      this.a0 = oneMinusAlpha * this.a0 + alpha * this.a;
      this.t0 = t;
    }
  }
};

const b2Vec3 = class {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  SetZero() {
    this.x = this.y = this.z = 0.0;
  }
  Set(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  SetV(v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
  }
  GetNegative() {
    return new b2Vec3(-this.x, -this.y, -this.z);
  }
  NegativeSelf() {
    this.x = -this.x;
    this.y = -this.y;
    this.z = -this.z;
  }
  Copy() {
    return new b2Vec3(this.x, this.y, this.z);
  }
  Add(v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
  }
  Subtract(v) {
    this.x -= v.x;
    this.y -= v.y;
    this.z -= v.z;
  }
  Multiply(a = 0) {
    this.x *= a;
    this.y *= a;
    this.z *= a;
  }
};

const b2AABB = class {
  constructor() {
    this.lowerBound = new b2Vec2();
    this.upperBound = new b2Vec2();
  }
  IsValid() {
    const dX = this.upperBound.x - this.lowerBound.x;
    const dY = this.upperBound.y - this.lowerBound.y;
    let valid = dX >= 0.0 && dY >= 0.0;
    valid = valid && this.lowerBound.IsValid() && this.upperBound.IsValid();
    return valid;
  }
  GetCenter() {
    return new b2Vec2((this.lowerBound.x + this.upperBound.x) / 2, (this.lowerBound.y + this.upperBound.y) / 2);
  }
  GetExtents() {
    return new b2Vec2((this.upperBound.x - this.lowerBound.x) / 2, (this.upperBound.y - this.lowerBound.y) / 2);
  }
  Contains(aabb) {
    return this.lowerBound.x <= aabb.lowerBound.x && this.lowerBound.y <= aabb.lowerBound.y && aabb.upperBound.x <= this.upperBound.x && aabb.upperBound.y <= this.upperBound.y;
  }
  RayCast(output, input) {
    let tmin = -Number.MAX_VALUE;
    let tmax = Number.MAX_VALUE;
    const pX = input.p1.x;
    const pY = input.p1.y;
    const dX = input.p2.x - input.p1.x;
    const dY = input.p2.y - input.p1.y;
    const absDX = Math.abs(dX);
    const absDY = Math.abs(dY);
    const normal = output.normal;
    let inv_d;
    let t1;
    let t2;
    let t3;
    let s;
    if (absDX < b2Settings.b2_epsilon) {
      if (pX < this.lowerBound.x || this.upperBound.x < pX) return false;
    } else {
      inv_d = 1.0 / dX;
      t1 = (this.lowerBound.x - pX) * inv_d;
      t2 = (this.upperBound.x - pX) * inv_d;
      s = -1.0;
      if (t1 > t2) {
        t3 = t1;
        t1 = t2;
        t2 = t3;
        s = 1.0;
      }
      if (t1 > tmin) {
        normal.x = s;
        normal.y = 0;
        tmin = t1;
      }
      tmax = Math.min(tmax, t2);
      if (tmin > tmax) return false;
    }
    if (absDY < b2Settings.b2_epsilon) {
      if (pY < this.lowerBound.y || this.upperBound.y < pY) return false;
    } else {
      inv_d = 1.0 / dY;
      t1 = (this.lowerBound.y - pY) * inv_d;
      t2 = (this.upperBound.y - pY) * inv_d;
      s = -1.0;
      if (t1 > t2) {
        t3 = t1;
        t1 = t2;
        t2 = t3;
        s = 1.0;
      }
      if (t1 > tmin) {
        normal.y = s;
        normal.x = 0;
        tmin = t1;
      }
      tmax = Math.min(tmax, t2);
      if (tmin > tmax) return false;
    }
    if (tmin < 0.0 || input.maxFraction < tmin) return false;
    output.fraction = tmin;
    return true;
  }
  TestOverlap(other) {
    const d1X = other.lowerBound.x - this.upperBound.x;
    const d1Y = other.lowerBound.y - this.upperBound.y;
    const d2X = this.lowerBound.x - other.upperBound.x;
    const d2Y = this.lowerBound.y - other.upperBound.y;
    if (d1X > 0.0 || d1Y > 0.0) return false;
    if (d2X > 0.0 || d2Y > 0.0) return false;
    return true;
  }
  static Combine(aabb1, aabb2) {
    const aabb = new b2AABB();
    aabb.Combine(aabb1, aabb2);
    return aabb;
  }
  Combine(aabb1, aabb2) {
    this.lowerBound.x = Math.min(aabb1.lowerBound.x, aabb2.lowerBound.x);
    this.lowerBound.y = Math.min(aabb1.lowerBound.y, aabb2.lowerBound.y);
    this.upperBound.x = Math.max(aabb1.upperBound.x, aabb2.upperBound.x);
    this.upperBound.y = Math.max(aabb1.upperBound.y, aabb2.upperBound.y);
  }
};

const b2Bound = class {
  IsLower() {
    return (this.value & 1) == 0;
  }
  IsUpper() {
    return (this.value & 1) == 1;
  }
  Swap(b) {
    const tempValue = this.value;
    const tempProxy = this.proxy;
    const tempStabbingCount = this.stabbingCount;
    this.value = b.value;
    this.proxy = b.proxy;
    this.stabbingCount = b.stabbingCount;
    b.value = tempValue;
    b.proxy = tempProxy;
    b.stabbingCount = tempStabbingCount;
  }
};

const b2BoundValues = class {
  constructor() {
    this.lowerValues = [];
    this.lowerValues[0] = 0.0;
    this.lowerValues[1] = 0.0;
    this.upperValues = [];
    this.upperValues[0] = 0.0;
    this.upperValues[1] = 0.0;
  }
};

const b2Manifold = class {
  static e_circles = 0x0001;
  static e_faceA = 0x0002;
  static e_faceB = 0x0004;
  constructor() {
    this.m_pointCount = 0;
    this.m_points = new Array(b2Settings.b2_maxManifoldPoints);
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.m_points[i] = new b2ManifoldPoint();
    }
    this.m_localPlaneNormal = new b2Vec2();
    this.m_localPoint = new b2Vec2();
  }
  Reset() {
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.m_points[i].Reset();
    }
    this.m_localPlaneNormal.SetZero();
    this.m_localPoint.SetZero();
    this.m_type = 0;
    this.m_pointCount = 0;
  }
  Set(m) {
    this.m_pointCount = m.m_pointCount;
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.m_points[i].Set(m.m_points[i]);
    }
    this.m_localPlaneNormal.SetV(m.m_localPlaneNormal);
    this.m_localPoint.SetV(m.m_localPoint);
    this.m_type = m.m_type;
  }
  Copy() {
    const copy = new b2Manifold();
    copy.Set(this);
    return copy;
  }
};

const b2ManifoldPoint = class {
  constructor() {
    this.m_localPoint = new b2Vec2();
    this.m_id = new b2ContactID();
    this.Reset();
  }
  Reset() {
    this.m_localPoint.SetZero();
    this.m_normalImpulse = 0.0;
    this.m_tangentImpulse = 0.0;
    this.m_id.key = 0;
  }
  Set(m) {
    this.m_localPoint.SetV(m.m_localPoint);
    this.m_normalImpulse = m.m_normalImpulse;
    this.m_tangentImpulse = m.m_tangentImpulse;
    this.m_id.Set(m.m_id);
  }
};

const b2ContactID = class {
  constructor() {
    this.features = new Features();
    this.features._m_id = this;
  }
  Set(id) {
    this.key = id._key;
  }
  Copy() {
    const id = new b2ContactID();
    id.key = this.key;
    return id;
  }
  get key() {
    return this._key;
  }
  set key(value = 0) {
    this._key = value;
    this.features._referenceEdge = this._key & 0x000000ff;
    this.features._incidentEdge = ((this._key & 0x0000ff00) >> 8) & 0x000000ff;
    this.features._incidentVertex = ((this._key & 0x00ff0000) >> 16) & 0x000000ff;
    this.features._flip = ((this._key & 0xff000000) >> 24) & 0x000000ff;
  }
};

const Features = class {
  get referenceEdge() {
    return this._referenceEdge;
  }
  set referenceEdge(value = 0) {
    this._referenceEdge = value;
    this._m_id._key = (this._m_id._key & 0xffffff00) | (this._referenceEdge & 0x000000ff);
  }
  get incidentEdge() {
    return this._incidentEdge;
  }
  set incidentEdge(value = 0) {
    this._incidentEdge = value;
    this._m_id._key = (this._m_id._key & 0xffff00ff) | ((this._incidentEdge << 8) & 0x0000ff00);
  }
  get incidentVertex() {
    return this._incidentVertex;
  }
  set incidentVertex(value = 0) {
    this._incidentVertex = value;
    this._m_id._key = (this._m_id._key & 0xff00ffff) | ((this._incidentVertex << 16) & 0x00ff0000);
  }
  get flip() {
    return this._flip;
  }
  set flip(value = 0) {
    this._flip = value;
    this._m_id._key = (this._m_id._key & 0x00ffffff) | ((this._flip << 24) & 0xff000000);
  }
};

const ClipVertex = class {
  constructor() {
    this.v = new b2Vec2();
    this.id = new b2ContactID();
  }
  Set(other) {
    this.v.SetV(other.v);
    this.id.Set(other.id);
  }
};

const b2Collision = class {
  static s_edgeAO = new Array(1).fill(0);
  static s_edgeBO = new Array(1).fill(0);
  static b2_nullFeature = 0x000000ff;
  static s_incidentEdge = this.MakeClipPointVector();
  static s_clipPoints1 = this.MakeClipPointVector();
  static s_clipPoints2 = this.MakeClipPointVector();
  static s_tangent = new b2Vec2();
  static s_tangent2 = new b2Vec2();
  static b2CollidePolyTempVec = new b2Vec2();
  static ClipSegmentToLine(vOut, vIn, normal, offset = 0) {
    let numOut = 0;
    let cv = vIn[0];
    const vIn0 = cv.v;
    cv = vIn[1];
    const vIn1 = cv.v;
    const distance0 = normal.x * vIn0.x + normal.y * vIn0.y - offset;
    const distance1 = normal.x * vIn1.x + normal.y * vIn1.y - offset;
    if (distance0 <= 0.0) vOut[numOut++].Set(vIn[0]);
    if (distance1 <= 0.0) vOut[numOut++].Set(vIn[1]);
    if (distance0 * distance1 < 0.0) {
      const interp = distance0 / (distance0 - distance1);
      cv = vOut[numOut];
      const tVec = cv.v;
      tVec.x = vIn0.x + interp * (vIn1.x - vIn0.x);
      tVec.y = vIn0.y + interp * (vIn1.y - vIn0.y);
      if (distance0 > 0.0) {
        cv.id = vIn[0].id;
      } else {
        cv.id = vIn[1].id;
      }
      ++numOut;
    }
    return numOut;
  }
  static EdgeSeparation(poly1, xf1, edge1 = 0, poly2, xf2) {
    const count2 = poly2.m_vertexCount;
    const vertices2 = poly2.m_vertices;
    let tMat = xf1.R;
    let tVec = poly1.m_normals[edge1];
    const normal1WorldX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    const normal1WorldY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    tMat = xf2.R;
    const normal1X = tMat.col1.x * normal1WorldX + tMat.col1.y * normal1WorldY;
    const normal1Y = tMat.col2.x * normal1WorldX + tMat.col2.y * normal1WorldY;
    let index = 0;
    let minDot = Number.MAX_VALUE;
    for (let i = 0; i < count2; ++i) {
      tVec = vertices2[i];
      const dot = tVec.x * normal1X + tVec.y * normal1Y;
      if (dot < minDot) {
        minDot = dot;
        index = i;
      }
    }
    tVec = vertices2[index];
    let v2X = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    let v2Y = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    tVec = poly1.m_vertices[edge1];
    tMat = xf1.R;
    v2X -= xf1.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    v2Y -= xf1.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    return v2X * normal1WorldX + v2Y * normal1WorldY;
  }
  static FindMaxSeparation(edgeIndex, poly1, xf1, poly2, xf2) {
    const count1 = poly1.m_vertexCount;
    const normals1 = poly1.m_normals;
    let tMat = xf2.R;
    let tVec = poly2.m_centroid;
    let dX = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    let dY = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    tMat = xf1.R;
    tVec = poly1.m_centroid;
    dX -= xf1.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    dY -= xf1.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    const dLocal1X = dX * tMat.col1.x + dY * tMat.col1.y;
    const dLocal1Y = dX * tMat.col2.x + dY * tMat.col2.y;
    let edge = 0;
    let maxDot = -Number.MAX_VALUE;
    for (let i = 0; i < count1; ++i) {
      tVec = normals1[i];
      const dot = tVec.x * dLocal1X + tVec.y * dLocal1Y;
      if (dot > maxDot) {
        maxDot = dot;
        edge = i;
      }
    }
    let s = b2Collision.EdgeSeparation(poly1, xf1, edge, poly2, xf2);
    const prevEdge = edge - 1 >= 0 ? edge - 1 : count1 - 1;
    const sPrev = b2Collision.EdgeSeparation(poly1, xf1, prevEdge, poly2, xf2);
    const nextEdge = edge + 1 < count1 ? edge + 1 : 0;
    const sNext = b2Collision.EdgeSeparation(poly1, xf1, nextEdge, poly2, xf2);
    let bestEdge;
    let bestSeparation;
    let increment;
    if (sPrev > s && sPrev > sNext) {
      increment = -1;
      bestEdge = prevEdge;
      bestSeparation = sPrev;
    } else if (sNext > s) {
      increment = 1;
      bestEdge = nextEdge;
      bestSeparation = sNext;
    } else {
      edgeIndex[0] = edge;
      return s;
    }
    while (true) {
      if (increment == -1) edge = bestEdge - 1 >= 0 ? bestEdge - 1 : count1 - 1;
      else edge = bestEdge + 1 < count1 ? bestEdge + 1 : 0;
      s = b2Collision.EdgeSeparation(poly1, xf1, edge, poly2, xf2);
      if (s > bestSeparation) {
        bestEdge = edge;
        bestSeparation = s;
      } else {
        break;
      }
    }
    edgeIndex[0] = bestEdge;
    return bestSeparation;
  }
  static FindIncidentEdge(c, poly1, xf1, edge1 = 0, poly2, xf2) {
    const normals1 = poly1.m_normals;
    const count2 = poly2.m_vertexCount;
    const vertices2 = poly2.m_vertices;
    const normals2 = poly2.m_normals;
    let tMat = xf1.R;
    let tVec = normals1[edge1];
    let normal1X = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    let normal1Y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    tMat = xf2.R;
    const tX = tMat.col1.x * normal1X + tMat.col1.y * normal1Y;
    normal1Y = tMat.col2.x * normal1X + tMat.col2.y * normal1Y;
    normal1X = tX;
    let index = 0;
    let minDot = Number.MAX_VALUE;
    for (let i = 0; i < count2; ++i) {
      tVec = normals2[i];
      const dot = normal1X * tVec.x + normal1Y * tVec.y;
      if (dot < minDot) {
        minDot = dot;
        index = i;
      }
    }
    const i1 = index;
    const i2 = i1 + 1 < count2 ? i1 + 1 : 0;
    let tClip = c[0];
    tVec = vertices2[i1];
    tMat = xf2.R;
    tClip.v.x = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    tClip.v.y = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    tClip.id.features.referenceEdge = edge1;
    tClip.id.features.incidentEdge = i1;
    tClip.id.features.incidentVertex = 0;
    tClip = c[1];
    tVec = vertices2[i2];
    tMat = xf2.R;
    tClip.v.x = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    tClip.v.y = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    tClip.id.features.referenceEdge = edge1;
    tClip.id.features.incidentEdge = i2;
    tClip.id.features.incidentVertex = 1;
  }
  static MakeClipPointVector() {
    const r = new Array(2);
    r[0] = new ClipVertex();
    r[1] = new ClipVertex();
    return r;
  }
  static CollidePolygons(manifold, polyA, xfA, polyB, xfB) {
    let cv;
    manifold.m_pointCount = 0;
    const totalRadius = polyA.m_radius + polyB.m_radius;
    let edgeA = 0;
    b2Collision.s_edgeAO[0] = edgeA;
    const separationA = b2Collision.FindMaxSeparation(b2Collision.s_edgeAO, polyA, xfA, polyB, xfB);
    edgeA = b2Collision.s_edgeAO[0];
    if (separationA > totalRadius) return;
    let edgeB = 0;
    b2Collision.s_edgeBO[0] = edgeB;
    const separationB = b2Collision.FindMaxSeparation(b2Collision.s_edgeBO, polyB, xfB, polyA, xfA);
    edgeB = b2Collision.s_edgeBO[0];
    if (separationB > totalRadius) return;
    let poly1;
    let poly2;
    let xf1;
    let xf2;
    let edge1;
    let flip;
    const k_relativeTol = 0.98;
    const k_absoluteTol = 0.001;
    if (separationB > k_relativeTol * separationA + k_absoluteTol) {
      poly1 = polyB;
      poly2 = polyA;
      xf1 = xfB;
      xf2 = xfA;
      edge1 = edgeB;
      manifold.m_type = b2Manifold.e_faceB;
      flip = 1;
    } else {
      poly1 = polyA;
      poly2 = polyB;
      xf1 = xfA;
      xf2 = xfB;
      edge1 = edgeA;
      manifold.m_type = b2Manifold.e_faceA;
      flip = 0;
    }
    const incidentEdge = b2Collision.s_incidentEdge;
    b2Collision.FindIncidentEdge(incidentEdge, poly1, xf1, edge1, poly2, xf2);
    const count1 = poly1.m_vertexCount;
    const vertices1 = poly1.m_vertices;
    const local_v11 = vertices1[edge1];
    let local_v12;
    if (edge1 + 1 < count1) {
      local_v12 = vertices1[edge1 + 1];
    } else {
      local_v12 = vertices1[0];
    }
    let localTangentX = local_v12.x - local_v11.x;
    let localTangentY = local_v12.y - local_v11.y;
    const localTangentLength = Math.sqrt(localTangentX * localTangentX + localTangentY * localTangentY);
    if (!(localTangentLength < b2Settings.b2_epsilon)) {
      const invLength = 1.0 / localTangentLength;
      localTangentX *= invLength;
      localTangentY *= invLength;
    }
    let tMat = xf1.R;
    const tangentX = tMat.col1.x * localTangentX + tMat.col2.x * localTangentY;
    const tangentY = tMat.col1.y * localTangentX + tMat.col2.y * localTangentY;
    const tangent = b2Collision.s_tangent;
    tangent.x = tangentX;
    tangent.y = tangentY;
    const tangent2 = b2Collision.s_tangent2;
    tangent2.x = -tangentX;
    tangent2.y = -tangentY;
    const normalX = tangentY;
    const normalY = -tangentX;
    const v11X = xf1.position.x + (tMat.col1.x * local_v11.x + tMat.col2.x * local_v11.y);
    const v11Y = xf1.position.y + (tMat.col1.y * local_v11.x + tMat.col2.y * local_v11.y);
    const v12X = xf1.position.x + (tMat.col1.x * local_v12.x + tMat.col2.x * local_v12.y);
    const v12Y = xf1.position.y + (tMat.col1.y * local_v12.x + tMat.col2.y * local_v12.y);
    const frontOffset = normalX * v11X + normalY * v11Y;
    const sideOffset1 = (-tangentX * v11X) - tangentY * v11Y + totalRadius;
    const sideOffset2 = tangentX * v12X + tangentY * v12Y + totalRadius;
    const clipPoints1 = b2Collision.s_clipPoints1;
    const clipPoints2 = b2Collision.s_clipPoints2;
    let np = b2Collision.ClipSegmentToLine(clipPoints1, incidentEdge, tangent2, sideOffset1);
    if (np < 2) return;
    np = b2Collision.ClipSegmentToLine(clipPoints2, clipPoints1, tangent, sideOffset2);
    if (np < 2) return;
    manifold.m_localPlaneNormal.x = localTangentY;
    manifold.m_localPlaneNormal.y = -localTangentX;
    manifold.m_localPoint.x = 0.5 * (local_v11.x + local_v12.x);
    manifold.m_localPoint.y = 0.5 * (local_v11.y + local_v12.y);
    tMat = xf2.R;
    let pointCount = 0;
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; ++i) {
      cv = clipPoints2[i];
      const separation = normalX * cv.v.x + normalY * cv.v.y - frontOffset;
      if (separation <= totalRadius) {
        const cp = manifold.m_points[pointCount];
        const tX = cv.v.x - xf2.position.x;
        const tY = cv.v.y - xf2.position.y;
        cp.m_localPoint.x = tX * tMat.col1.x + tY * tMat.col1.y;
        cp.m_localPoint.y = tX * tMat.col2.x + tY * tMat.col2.y;
        cp.m_id.Set(cv.id);
        cp.m_id.features.flip = flip;
        ++pointCount;
      }
    }
    manifold.m_pointCount = pointCount;
  }
  static CollideCircles(manifold, circle1, xf1, circle2, xf2) {
    manifold.m_pointCount = 0;
    let tMat = xf1.R;
    let tVec = circle1.m_p;
    const p1X = xf1.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    const p1Y = xf1.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    tMat = xf2.R;
    tVec = circle2.m_p;
    const p2X = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    const p2Y = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    const distSqr = dX * dX + dY * dY;
    const radius = circle1.m_radius + circle2.m_radius;
    if (distSqr > radius * radius) {
      return;
    }
    manifold.m_type = b2Manifold.e_circles;
    manifold.m_localPoint.SetV(circle1.m_p);
    manifold.m_localPlaneNormal.SetZero();
    manifold.m_pointCount = 1;
    manifold.m_points[0].m_localPoint.SetV(circle2.m_p);
    manifold.m_points[0].m_id.key = 0;
  }
  static CollidePolygonAndCircle(manifold, polygon, xf1, circle, xf2) {
    manifold.m_pointCount = 0;
    let tMat = xf2.R;
    let tVec = circle.m_p;
    const cX = xf2.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    const cY = xf2.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    let dX = cX - xf1.position.x;
    let dY = cY - xf1.position.y;
    tMat = xf1.R;
    const cLocalX = dX * tMat.col1.x + dY * tMat.col1.y;
    const cLocalY = dX * tMat.col2.x + dY * tMat.col2.y;
    let normalIndex = 0;
    let separation = -Number.MAX_VALUE;
    const radius = polygon.m_radius + circle.m_radius;
    const vertexCount = polygon.m_vertexCount;
    const vertices = polygon.m_vertices;
    const normals = polygon.m_normals;
    for (let i = 0; i < vertexCount; ++i) {
      tVec = vertices[i];
      dX = cLocalX - tVec.x;
      dY = cLocalY - tVec.y;
      tVec = normals[i];
      const s = tVec.x * dX + tVec.y * dY;
      if (s > radius) {
        return;
      }
      if (s > separation) {
        separation = s;
        normalIndex = i;
      }
    }
    const vertIndex1 = normalIndex;
    const vertIndex2 = vertIndex1 + 1 < vertexCount ? vertIndex1 + 1 : 0;
    const v1 = vertices[vertIndex1];
    const v2 = vertices[vertIndex2];
    if (separation < b2Settings.b2_epsilon) {
      manifold.m_pointCount = 1;
      manifold.m_type = b2Manifold.e_faceA;
      manifold.m_localPlaneNormal.SetV(normals[vertIndex1]);
      manifold.m_localPoint.x = 0.5 * (v1.x + v2.x);
      manifold.m_localPoint.y = 0.5 * (v1.y + v2.y);
      manifold.m_points[0].m_localPoint.SetV(circle.m_p);
      manifold.m_points[0].m_id.key = 0;
      return;
    }
    const u1 = (cLocalX - v1.x) * (v2.x - v1.x) + (cLocalY - v1.y) * (v2.y - v1.y);
    const u2 = (cLocalX - v2.x) * (v1.x - v2.x) + (cLocalY - v2.y) * (v1.y - v2.y);
    if (u1 <= 0.0) {
      if ((cLocalX - v1.x) * (cLocalX - v1.x) + (cLocalY - v1.y) * (cLocalY - v1.y) > radius * radius) return;
      manifold.m_pointCount = 1;
      manifold.m_type = b2Manifold.e_faceA;
      manifold.m_localPlaneNormal.x = cLocalX - v1.x;
      manifold.m_localPlaneNormal.y = cLocalY - v1.y;
      manifold.m_localPlaneNormal.Normalize();
      manifold.m_localPoint.SetV(v1);
      manifold.m_points[0].m_localPoint.SetV(circle.m_p);
      manifold.m_points[0].m_id.key = 0;
    } else if (u2 <= 0) {
      if ((cLocalX - v2.x) * (cLocalX - v2.x) + (cLocalY - v2.y) * (cLocalY - v2.y) > radius * radius) return;
      manifold.m_pointCount = 1;
      manifold.m_type = b2Manifold.e_faceA;
      manifold.m_localPlaneNormal.x = cLocalX - v2.x;
      manifold.m_localPlaneNormal.y = cLocalY - v2.y;
      manifold.m_localPlaneNormal.Normalize();
      manifold.m_localPoint.SetV(v2);
      manifold.m_points[0].m_localPoint.SetV(circle.m_p);
      manifold.m_points[0].m_id.key = 0;
    } else {
      const faceCenterX = 0.5 * (v1.x + v2.x);
      const faceCenterY = 0.5 * (v1.y + v2.y);
      tVec = normals[vertIndex1];
      separation = (cLocalX - faceCenterX) * tVec.x + (cLocalY - faceCenterY) * tVec.y;
      if (separation > radius) return;
      manifold.m_pointCount = 1;
      manifold.m_type = b2Manifold.e_faceA;
      manifold.m_localPlaneNormal.x = tVec.x;
      manifold.m_localPlaneNormal.y = tVec.y;
      manifold.m_localPlaneNormal.Normalize();
      manifold.m_localPoint.Set(faceCenterX, faceCenterY);
      manifold.m_points[0].m_localPoint.SetV(circle.m_p);
      manifold.m_points[0].m_id.key = 0;
    }
  }
  static TestOverlap(a, b) {
    let t1 = b.lowerBound;
    let t2 = a.upperBound;
    const d1X = t1.x - t2.x;
    const d1Y = t1.y - t2.y;
    t1 = a.lowerBound;
    t2 = b.upperBound;
    const d2X = t1.x - t2.x;
    const d2Y = t1.y - t2.y;
    if (d1X > 0.0 || d1Y > 0.0) return false;
    if (d2X > 0.0 || d2Y > 0.0) return false;
    return true;
  }
};

const b2ContactPoint = class {
  constructor() {
    this.position = new b2Vec2();
    this.velocity = new b2Vec2();
    this.normal = new b2Vec2();
    this.id = new b2ContactID();
  }
};

const b2Simplex = class {
  constructor() {
    this.m_v1 = new b2SimplexVertex();
    this.m_v2 = new b2SimplexVertex();
    this.m_v3 = new b2SimplexVertex();
    this.m_vertices = new Array(3);
    this.m_vertices[0] = this.m_v1;
    this.m_vertices[1] = this.m_v2;
    this.m_vertices[2] = this.m_v3;
  }
  ReadCache(cache, proxyA, transformA, proxyB, transformB) {
    const count = cache.count;
    b2Settings.b2Assert(0 <= count && count <= 3);
    let wALocal;
    let wBLocal;
    let tMat;
    this.m_count = count;
    const vertices = this.m_vertices;
    for (let i = 0; i < count; i++) {
      const v = vertices[i];
      v.indexA = cache.indexA[i];
      v.indexB = cache.indexB[i];
      wALocal = proxyA.GetVertex(v.indexA);
      wBLocal = proxyB.GetVertex(v.indexB);
      tMat = transformA.R;
      v.wA.x = tMat.col1.x * wALocal.x + tMat.col2.x * wALocal.y + transformA.position.x;
      v.wA.y = tMat.col1.y * wALocal.x + tMat.col2.y * wALocal.y + transformA.position.y;
      tMat = transformB.R;
      v.wB.x = tMat.col1.x * wBLocal.x + tMat.col2.x * wBLocal.y + transformB.position.x;
      v.wB.y = tMat.col1.y * wBLocal.x + tMat.col2.y * wBLocal.y + transformB.position.y;
      v.w.x = v.wB.x - v.wA.x;
      v.w.y = v.wB.y - v.wA.y;
      v.a = 0;
    }
    if (count > 1) {
      const metric1 = cache.metric;
      const metric2 = this.GetMetric();
      if (metric2 < 0.5 * metric1 || 2.0 * metric1 < metric2 || metric2 < b2Settings.b2_epsilon) {
        this.m_count = 0;
      }
    }
    if (this.m_count == 0) {
      const v = vertices[0];
      v.indexA = 0;
      v.indexB = 0;
      wALocal = proxyA.GetVertex(0);
      wBLocal = proxyB.GetVertex(0);
      tMat = transformA.R;
      v.wA.x = tMat.col1.x * wALocal.x + tMat.col2.x * wALocal.y + transformA.position.x;
      v.wA.y = tMat.col1.y * wALocal.x + tMat.col2.y * wALocal.y + transformA.position.y;
      tMat = transformB.R;
      v.wB.x = tMat.col1.x * wBLocal.x + tMat.col2.x * wBLocal.y + transformB.position.x;
      v.wB.y = tMat.col1.y * wBLocal.x + tMat.col2.y * wBLocal.y + transformB.position.y;
      v.w.x = v.wB.x - v.wA.x;
      v.w.y = v.wB.y - v.wA.y;
      this.m_count = 1;
    }
  }
  WriteCache(cache) {
    cache.metric = this.GetMetric();
    cache.count = this.m_count;
    const vertices = this.m_vertices;
    for (let i = 0; i < this.m_count; i++) {
      cache.indexA[i] = vertices[i].indexA;
      cache.indexB[i] = vertices[i].indexB;
    }
  }
  GetSearchDirection() {
    switch (this.m_count) {
    case 1:
      return this.m_v1.w.GetNegative();
    case 2: {
      const e12 = b2Math.SubtractVV(this.m_v2.w, this.m_v1.w);
      const sgn = b2Math.CrossVV(e12, this.m_v1.w.GetNegative());
      if (sgn > 0.0) {
        return b2Math.CrossFV(1.0, e12);
      } else {
        return b2Math.CrossVF(e12, 1.0);
      }
    }
    default:
      b2Settings.b2Assert(false);
      return new b2Vec2();
    }
  }
  GetClosestPoint() {
    switch (this.m_count) {
    case 0:
      b2Settings.b2Assert(false);
      return new b2Vec2();
    case 1:
      return this.m_v1.w;
    case 2:
      return new b2Vec2(this.m_v1.a * this.m_v1.w.x + this.m_v2.a * this.m_v2.w.x, this.m_v1.a * this.m_v1.w.y + this.m_v2.a * this.m_v2.w.y);
    default:
      b2Settings.b2Assert(false);
      return new b2Vec2();
    }
  }
  GetWitnessPoints(pA, pB) {
    switch (this.m_count) {
    case 0:
      b2Settings.b2Assert(false);
      break;
    case 1:
      pA.SetV(this.m_v1.wA);
      pB.SetV(this.m_v1.wB);
      break;
    case 2:
      pA.x = this.m_v1.a * this.m_v1.wA.x + this.m_v2.a * this.m_v2.wA.x;
      pA.y = this.m_v1.a * this.m_v1.wA.y + this.m_v2.a * this.m_v2.wA.y;
      pB.x = this.m_v1.a * this.m_v1.wB.x + this.m_v2.a * this.m_v2.wB.x;
      pB.y = this.m_v1.a * this.m_v1.wB.y + this.m_v2.a * this.m_v2.wB.y;
      break;
    case 3:
      pB.x = pA.x = this.m_v1.a * this.m_v1.wA.x + this.m_v2.a * this.m_v2.wA.x + this.m_v3.a * this.m_v3.wA.x;
      pB.y = pA.y = this.m_v1.a * this.m_v1.wA.y + this.m_v2.a * this.m_v2.wA.y + this.m_v3.a * this.m_v3.wA.y;
      break;
    default:
      b2Settings.b2Assert(false);
      break;
    }
  }
  GetMetric() {
    switch (this.m_count) {
    case 0:
      b2Settings.b2Assert(false);
      return 0.0;
    case 1:
      return 0.0;
    case 2: {
      const dX = this.m_v1.w.x - this.m_v2.w.x;
      const dY = this.m_v1.w.y - this.m_v2.w.y;
      return Math.sqrt(dX * dX + dY * dY);
    }
    case 3: {
      const e12X = this.m_v2.w.x - this.m_v1.w.x;
      const e12Y = this.m_v2.w.y - this.m_v1.w.y;
      const e13X = this.m_v3.w.x - this.m_v1.w.x;
      const e13Y = this.m_v3.w.y - this.m_v1.w.y;
      return e12X * e13Y - e12Y * e13X;
    }
    default:
      b2Settings.b2Assert(false);
      return 0.0;
    }
  }
  Solve2() {
    const w1 = this.m_v1.w;
    const w2 = this.m_v2.w;
    const e12X = w2.x - w1.x;
    const e12Y = w2.y - w1.y;
    const d12_2 = -(w1.x * e12X + w1.y * e12Y);
    if (d12_2 <= 0.0) {
      this.m_v1.a = 1.0;
      this.m_count = 1;
      return;
    }
    const d12_1 = w2.x * e12X + w2.y * e12Y;
    if (d12_1 <= 0.0) {
      this.m_v2.a = 1.0;
      this.m_count = 1;
      this.m_v1.Set(this.m_v2);
      return;
    }
    const inv_d12 = 1.0 / (d12_1 + d12_2);
    this.m_v1.a = d12_1 * inv_d12;
    this.m_v2.a = d12_2 * inv_d12;
    this.m_count = 2;
  }
  Solve3() {
    const w1 = this.m_v1.w;
    const w2 = this.m_v2.w;
    const w3 = this.m_v3.w;
    const e12X = w2.x - w1.x;
    const e12Y = w2.y - w1.y;
    const w1e12 = w1.x * e12X + w1.y * e12Y;
    const w2e12 = w2.x * e12X + w2.y * e12Y;
    const d12_1 = w2e12;
    const d12_2 = -w1e12;
    const e13X = w3.x - w1.x;
    const e13Y = w3.y - w1.y;
    const w1e13 = w1.x * e13X + w1.y * e13Y;
    const w3e13 = w3.x * e13X + w3.y * e13Y;
    const d13_1 = w3e13;
    const d13_2 = -w1e13;
    const e23X = w3.x - w2.x;
    const e23Y = w3.y - w2.y;
    const w2e23 = w2.x * e23X + w2.y * e23Y;
    const w3e23 = w3.x * e23X + w3.y * e23Y;
    const d23_1 = w3e23;
    const d23_2 = -w2e23;
    const n123 = e12X * e13Y - e12Y * e13X;
    const d123_1 = n123 * (w2.x * w3.y - w2.y * w3.x);
    const d123_2 = n123 * (w3.x * w1.y - w3.y * w1.x);
    const d123_3 = n123 * (w1.x * w2.y - w1.y * w2.x);
    if (d12_2 <= 0.0 && d13_2 <= 0.0) {
      this.m_v1.a = 1.0;
      this.m_count = 1;
      return;
    }
    if (d12_1 > 0.0 && d12_2 > 0.0 && d123_3 <= 0.0) {
      const inv_d12 = 1.0 / (d12_1 + d12_2);
      this.m_v1.a = d12_1 * inv_d12;
      this.m_v2.a = d12_2 * inv_d12;
      this.m_count = 2;
      return;
    }
    if (d13_1 > 0.0 && d13_2 > 0.0 && d123_2 <= 0.0) {
      const inv_d13 = 1.0 / (d13_1 + d13_2);
      this.m_v1.a = d13_1 * inv_d13;
      this.m_v3.a = d13_2 * inv_d13;
      this.m_count = 2;
      this.m_v2.Set(this.m_v3);
      return;
    }
    if (d12_1 <= 0.0 && d23_2 <= 0.0) {
      this.m_v2.a = 1.0;
      this.m_count = 1;
      this.m_v1.Set(this.m_v2);
      return;
    }
    if (d13_1 <= 0.0 && d23_1 <= 0.0) {
      this.m_v3.a = 1.0;
      this.m_count = 1;
      this.m_v1.Set(this.m_v3);
      return;
    }
    if (d23_1 > 0.0 && d23_2 > 0.0 && d123_1 <= 0.0) {
      const inv_d23 = 1.0 / (d23_1 + d23_2);
      this.m_v2.a = d23_1 * inv_d23;
      this.m_v3.a = d23_2 * inv_d23;
      this.m_count = 2;
      this.m_v1.Set(this.m_v3);
      return;
    }
    const inv_d123 = 1.0 / (d123_1 + d123_2 + d123_3);
    this.m_v1.a = d123_1 * inv_d123;
    this.m_v2.a = d123_2 * inv_d123;
    this.m_v3.a = d123_3 * inv_d123;
    this.m_count = 3;
  }
};

const b2SimplexVertex = class {
  constructor() {
    this.wA = new b2Vec2();
    this.wB = new b2Vec2();
    this.w = new b2Vec2();
  }
  Set(other) {
    this.wA.SetV(other.wA);
    this.wB.SetV(other.wB);
    this.w.SetV(other.w);
    this.a = other.a;
    this.indexA = other.indexA;
    this.indexB = other.indexB;
  }
};

const b2Distance = class {
  static s_saveA = new Array(3).fill(0);
  static s_saveB = new Array(3).fill(0);
  static s_simplex = new b2Simplex();
  static s_support = new b2Vec2();
  static s_normal = new b2Vec2();
  static Distance(output, cache, input) {
    ++b2Distance.b2_gjkCalls;
    const proxyA = input.proxyA;
    const proxyB = input.proxyB;
    const transformA = input.transformA;
    const transformB = input.transformB;
    const simplex = b2Distance.s_simplex;
    simplex.ReadCache(cache, proxyA, transformA, proxyB, transformB);
    const vertices = simplex.m_vertices;
    const w1 = simplex.m_v1.w;
    const w2 = simplex.m_v2.w;
    const support = b2Distance.s_support;
    const k_maxIters = 20;
    const saveA = b2Distance.s_saveA;
    const saveB = b2Distance.s_saveB;
    let saveCount;
    let tMat;
    let tVec;
    let iter = 0;
    while (iter < k_maxIters) {
      saveCount = simplex.m_count;
      for (let i = 0; i < saveCount; i++) {
        saveA[i] = vertices[i].indexA;
        saveB[i] = vertices[i].indexB;
      }
      switch (simplex.m_count) {
      case 1:
        break;
      case 2:
        simplex.Solve2();
        break;
      case 3:
        simplex.Solve3();
        break;
      default:
        b2Settings.b2Assert(false);
      }
      if (simplex.m_count == 3) {
        break;
      }
      let dX;
      let dY;
      if (simplex.m_count == 1) {
        dX = -w1.x;
        dY = -w1.y;
      } else {
        const e12X = w2.x - w1.x;
        const e12Y = w2.y - w1.y;
        if (e12X * -w1.y - e12Y * -w1.x > 0.0) {
          dX = -1.0 * e12Y;
          dY = 1.0 * e12X;
        } else {
          dX = 1.0 * e12Y;
          dY = -1.0 * e12X;
        }
      }
      if (dX * dX + dY * dY < b2Settings.b2_epsilon * b2Settings.b2_epsilon) {
        break;
      }
      const vertex = vertices[simplex.m_count];
      tMat = transformA.R;
      support.x = -dX * tMat.col1.x + -dY * tMat.col1.y;
      support.y = -dX * tMat.col2.x + -dY * tMat.col2.y;
      const indexA = proxyA.GetSupport(support);
      vertex.indexA = indexA;
      tVec = proxyA.GetVertex(indexA);
      const wAX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + transformA.position.x;
      const wAY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + transformA.position.y;
      tMat = transformB.R;
      support.x = dX * tMat.col1.x + dY * tMat.col1.y;
      support.y = dX * tMat.col2.x + dY * tMat.col2.y;
      const indexB = proxyB.GetSupport(support);
      vertex.indexB = indexB;
      tVec = proxyB.GetVertex(indexB);
      const wBX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + transformB.position.x;
      const wBY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + transformB.position.y;
      vertex.wA.x = wAX;
      vertex.wA.y = wAY;
      vertex.wB.x = wBX;
      vertex.wB.y = wBY;
      vertex.w.x = wBX - wAX;
      vertex.w.y = wBY - wAY;
      ++iter;
      ++b2Distance.b2_gjkIters;
      let duplicate = false;
      for (let i = 0; i < saveCount; i++) {
        if (indexA == saveA[i] && indexB == saveB[i]) {
          duplicate = true;
          break;
        }
      }
      if (duplicate) {
        break;
      }
      ++simplex.m_count;
    }
    b2Distance.b2_gjkMaxIters = b2Math.Max(b2Distance.b2_gjkMaxIters, iter);
    simplex.GetWitnessPoints(output.pointA, output.pointB);
    const distanceX = output.pointA.x - output.pointB.x;
    const distanceY = output.pointA.y - output.pointB.y;
    output.distance = Math.sqrt(distanceX * distanceX + distanceY * distanceY);
    output.iterations = iter;
    simplex.WriteCache(cache);
    if (input.useRadii) {
      const rA = proxyA.m_radius;
      const rB = proxyB.m_radius;
      if (output.distance > rA + rB && output.distance > b2Settings.b2_epsilon) {
        output.distance -= rA + rB;
        const normal = b2Distance.s_normal;
        normal.Set(output.pointB.x - output.pointA.x, output.pointB.y - output.pointA.y);
        normal.Normalize();
        output.pointA.x += rA * normal.x;
        output.pointA.y += rA * normal.y;
        output.pointB.x -= rB * normal.x;
        output.pointB.y -= rB * normal.y;
      } else {
        const pX = 0.5 * (output.pointA.x + output.pointB.x);
        const pY = 0.5 * (output.pointA.y + output.pointB.y);
        output.pointA.x = output.pointB.x = pX;
        output.pointA.y = output.pointB.y = pY;
        output.distance = 0.0;
      }
    }
  }
};

const b2DistanceInput = class {
};

const b2DistanceOutput = class {
  constructor() {
    this.pointA = new b2Vec2();
    this.pointB = new b2Vec2();
  }
};

const b2DistanceProxy = class {
  constructor() {
    this.m_circleVertices = [null];
  }
  Set(shape) {
    switch (shape.GetType()) {
    case b2Shape.e_circleShape: {
      const circle = shape;
      this.m_vertices = this.m_circleVertices;
      this.m_vertices[0] = circle.m_p;
      this.m_count = 1;
      this.m_radius = circle.m_radius;
    }
      break;
    case b2Shape.e_polygonShape: {
      const polygon = shape;
      this.m_vertices = polygon.m_vertices;
      this.m_count = polygon.m_vertexCount;
      this.m_radius = polygon.m_radius;
    }
      break;
    default:
      b2Settings.b2Assert(false);
    }
  }
  GetSupport(d) {
    let bestIndex = 0;
    let bestValue = this.m_vertices[0].x * d.x + this.m_vertices[0].y * d.y;
    for (let i = 1; i < this.m_count; ++i) {
      const value = this.m_vertices[i].x * d.x + this.m_vertices[i].y * d.y;
      if (value > bestValue) {
        bestIndex = i;
        bestValue = value;
      }
    }
    return bestIndex;
  }
  GetSupportVertex(d) {
    let bestIndex = 0;
    let bestValue = this.m_vertices[0].x * d.x + this.m_vertices[0].y * d.y;
    for (let i = 1; i < this.m_count; ++i) {
      const value = this.m_vertices[i].x * d.x + this.m_vertices[i].y * d.y;
      if (value > bestValue) {
        bestIndex = i;
        bestValue = value;
      }
    }
    return this.m_vertices[bestIndex];
  }
  GetVertexCount() {
    return this.m_count;
  }
  GetVertex(index = 0) {
    b2Settings.b2Assert(0 <= index && index < this.m_count);
    return this.m_vertices[index];
  }
};

const b2DynamicTree = class {
  constructor() {
    this.m_root = null;
    this.m_freeList = null;
    this.m_path = 0;
    this.m_insertionCount = 0;
    this.m_queryStacks = [];
  }
  CreateProxy(aabb, userData) {
    const node = this.AllocateNode();
    const extendX = b2Settings.b2_aabbExtension;
    const extendY = b2Settings.b2_aabbExtension;
    node.aabb.lowerBound.x = aabb.lowerBound.x - extendX;
    node.aabb.lowerBound.y = aabb.lowerBound.y - extendY;
    node.aabb.upperBound.x = aabb.upperBound.x + extendX;
    node.aabb.upperBound.y = aabb.upperBound.y + extendY;
    node.userData = userData;
    this.InsertLeaf(node);
    return node;
  }
  DestroyProxy(proxy) {
    this.RemoveLeaf(proxy);
    this.FreeNode(proxy);
  }
  MoveProxy(proxy, aabb, displacement) {
    b2Settings.b2Assert(proxy.IsLeaf());
    if (proxy.aabb.Contains(aabb)) {
      return false;
    }
    this.RemoveLeaf(proxy);
    const extendX = b2Settings.b2_aabbExtension + b2Settings.b2_aabbMultiplier * (displacement.x > 0 ? displacement.x : -displacement.x);
    const extendY = b2Settings.b2_aabbExtension + b2Settings.b2_aabbMultiplier * (displacement.y > 0 ? displacement.y : -displacement.y);
    proxy.aabb.lowerBound.x = aabb.lowerBound.x - extendX;
    proxy.aabb.lowerBound.y = aabb.lowerBound.y - extendY;
    proxy.aabb.upperBound.x = aabb.upperBound.x + extendX;
    proxy.aabb.upperBound.y = aabb.upperBound.y + extendY;
    this.InsertLeaf(proxy);
    return true;
  }
  Rebalance(iterations = 0) {
    if (this.m_root == null) return;
    for (let i = 0; i < iterations; i++) {
      let node = this.m_root;
      let bit = 0;
      while (node.IsLeaf() == false) {
        node = (this.m_path >> bit) & 1 ? node.child2 : node.child1;
        bit = (bit + 1) & 31;
      }
      ++this.m_path;
      this.RemoveLeaf(node);
      this.InsertLeaf(node);
    }
  }
  GetFatAABB(proxy) {
    return proxy.aabb;
  }
  GetUserData(proxy) {
    return proxy.userData;
  }
  Query(callback, aabb) {
    if (this.m_root == null) return;
    const stack = this.m_queryStacks.length > 0 ? this.m_queryStacks.pop() : [];
    let count = 0;
    stack[count++] = this.m_root;
    while (count > 0) {
      const node = stack[--count];
      if (node.aabb.TestOverlap(aabb)) {
        if (node.IsLeaf()) {
          const proceed = callback(node);
          if (!proceed) break;
        } else {
          stack[count++] = node.child1;
          stack[count++] = node.child2;
        }
      }
    }
    this.m_queryStacks.push(stack);
  }
  RayCast(callback, input) {
    if (this.m_root == null) return;
    const p1 = input.p1;
    const p2 = input.p2;
    let rX = p1.x - p2.x;
    let rY = p1.y - p2.y;
    const rLength = Math.sqrt(rX * rX + rY * rY);
    if (!(rLength < b2Settings.b2_epsilon)) {
      const invLength = 1.0 / rLength;
      rX *= invLength;
      rY *= invLength;
    }
    const vX = -1.0 * rY;
    const vY = 1.0 * rX;
    const absVX = b2Math.Abs(vX);
    const absVY = b2Math.Abs(vY);
    let maxFraction = input.maxFraction;
    let tX = p1.x + maxFraction * (p2.x - p1.x);
    let tY = p1.y + maxFraction * (p2.y - p1.y);
    let segmentLowerX = Math.min(p1.x, tX);
    let segmentLowerY = Math.min(p1.y, tY);
    let segmentUpperX = Math.max(p1.x, tX);
    let segmentUpperY = Math.max(p1.y, tY);
    let subInput = null;
    const stack = this.m_queryStacks.length > 0 ? this.m_queryStacks.pop() : [];
    let count = 0;
    stack[count++] = this.m_root;
    while (count > 0) {
      const node = stack[--count];
      const aabb = node.aabb;
      if (segmentLowerX - aabb.upperBound.x > 0.0 || segmentLowerY - aabb.upperBound.y > 0.0) continue;
      if (aabb.lowerBound.x - segmentUpperX > 0.0 || aabb.lowerBound.y - segmentUpperY > 0.0) continue;
      const cX = (aabb.lowerBound.x + aabb.upperBound.x) / 2;
      const cY = (aabb.lowerBound.y + aabb.upperBound.y) / 2;
      const hX = (aabb.upperBound.x - aabb.lowerBound.x) / 2;
      const hY = (aabb.upperBound.y - aabb.lowerBound.y) / 2;
      const separation = Math.abs(vX * (p1.x - cX) + vY * (p1.y - cY)) - absVX * hX - absVY * hY;
      if (separation > 0.0) continue;
      if (node.IsLeaf()) {
        if (subInput == null) subInput = new b2RayCastInput();
        subInput.p1 = input.p1;
        subInput.p2 = input.p2;
        subInput.maxFraction = maxFraction;
        maxFraction = callback(subInput, node);
        if (maxFraction == 0.0) break;
        if (maxFraction > 0.0) {
          tX = p1.x + maxFraction * (p2.x - p1.x);
          tY = p1.y + maxFraction * (p2.y - p1.y);
          segmentLowerX = Math.min(p1.x, tX);
          segmentLowerY = Math.min(p1.y, tY);
          segmentUpperX = Math.max(p1.x, tX);
          segmentUpperY = Math.max(p1.y, tY);
        }
      } else {
        stack[count++] = node.child1;
        stack[count++] = node.child2;
      }
    }
    this.m_queryStacks.push(stack);
  }
  AllocateNode() {
    if (this.m_freeList) {
      const node = this.m_freeList;
      this.m_freeList = node.parent;
      node.parent = null;
      node.child1 = null;
      node.child2 = null;
      return node;
    }
    return new b2DynamicTreeNode();
  }
  FreeNode(node) {
    node.parent = this.m_freeList;
    this.m_freeList = node;
  }
  InsertLeaf(leaf) {
    ++this.m_insertionCount;
    if (this.m_root == null) {
      this.m_root = leaf;
      this.m_root.parent = null;
      return;
    }
    const centerX = (leaf.aabb.lowerBound.x + leaf.aabb.upperBound.x) / 2;
    const centerY = (leaf.aabb.lowerBound.y + leaf.aabb.upperBound.y) / 2;
    let sibling = this.m_root;
    if (sibling.IsLeaf() == false) {
      do {
        const child1 = sibling.child1;
        const child2 = sibling.child2;
        const norm1 = Math.abs((child1.aabb.lowerBound.x + child1.aabb.upperBound.x) / 2 - centerX) + Math.abs((child1.aabb.lowerBound.y + child1.aabb.upperBound.y) / 2 - centerY);
        const norm2 = Math.abs((child2.aabb.lowerBound.x + child2.aabb.upperBound.x) / 2 - centerX) + Math.abs((child2.aabb.lowerBound.y + child2.aabb.upperBound.y) / 2 - centerY);
        if (norm1 < norm2) {
          sibling = child1;
        } else {
          sibling = child2;
        }
      }
      while (sibling.IsLeaf() == false)
    }
    let node1 = sibling.parent;
    let node2 = this.AllocateNode();
    node2.parent = node1;
    node2.userData = null;
    node2.aabb.Combine(leaf.aabb, sibling.aabb);
    if (node1) {
      if (sibling.parent.child1 == sibling) {
        node1.child1 = node2;
      } else {
        node1.child2 = node2;
      }
      node2.child1 = sibling;
      node2.child2 = leaf;
      sibling.parent = node2;
      leaf.parent = node2;
      do {
        if (node1.aabb.Contains(node2.aabb)) break;
        node1.aabb.Combine(node1.child1.aabb, node1.child2.aabb);
        node2 = node1;
        node1 = node1.parent;
      }
      while (node1)
    } else {
      node2.child1 = sibling;
      node2.child2 = leaf;
      sibling.parent = node2;
      leaf.parent = node2;
      this.m_root = node2;
    }
  }
  RemoveLeaf(leaf) {
    if (leaf == this.m_root) {
      this.m_root = null;
      return;
    }
    const node2 = leaf.parent;
    let node1 = node2.parent;
    let sibling;
    if (node2.child1 == leaf) {
      sibling = node2.child2;
    } else {
      sibling = node2.child1;
    }
    if (node1) {
      if (node1.child1 == node2) {
        node1.child1 = sibling;
      } else {
        node1.child2 = sibling;
      }
      sibling.parent = node1;
      this.FreeNode(node2);
      while (node1) {
        const aabb = node1.aabb;
        const lowerX = Math.min(node1.child1.aabb.lowerBound.x, node1.child2.aabb.lowerBound.x);
        const lowerY = Math.min(node1.child1.aabb.lowerBound.y, node1.child2.aabb.lowerBound.y);
        const upperX = Math.max(node1.child1.aabb.upperBound.x, node1.child2.aabb.upperBound.x);
        const upperY = Math.max(node1.child1.aabb.upperBound.y, node1.child2.aabb.upperBound.y);
        const contains = aabb.lowerBound.x <= lowerX && aabb.lowerBound.y <= lowerY && upperX <= aabb.upperBound.x && upperY <= aabb.upperBound.y;
        aabb.lowerBound.x = lowerX;
        aabb.lowerBound.y = lowerY;
        aabb.upperBound.x = upperX;
        aabb.upperBound.y = upperY;
        if (contains) break;
        node1 = node1.parent;
      }
    } else {
      this.m_root = sibling;
      sibling.parent = null;
      this.FreeNode(node2);
    }
  }
};

const b2DynamicTreeBroadPhase = class {
  constructor() {
    this.m_tree = new b2DynamicTree();
    this.m_moveBuffer = [];
    this.m_pairBuffer = [];
    this.m_pairCount = 0;
    this.m_proxyCount = 0;
    this.m_queryProxy = null;
    this.m_queryCallback = (proxy) => {
      if (proxy == this.m_queryProxy) return true;
      if (this.m_pairCount == this.m_pairBuffer.length) {
        this.m_pairBuffer[this.m_pairCount] = new b2DynamicTreePair();
      }
      const pair = this.m_pairBuffer[this.m_pairCount];
      pair.proxyA = this.m_queryProxy;
      pair.proxyB = proxy;
      ++this.m_pairCount;
      return true;
    };
  }
  CreateProxy(aabb, userData) {
    const proxy = this.m_tree.CreateProxy(aabb, userData);
    ++this.m_proxyCount;
    this.BufferMove(proxy);
    return proxy;
  }
  DestroyProxy(proxy) {
    this.UnBufferMove(proxy);
    --this.m_proxyCount;
    this.m_tree.DestroyProxy(proxy);
  }
  MoveProxy(proxy, aabb, displacement) {
    const buffer = this.m_tree.MoveProxy(proxy, aabb, displacement);
    if (buffer) {
      this.BufferMove(proxy);
    }
  }
  TestOverlap(proxyA, proxyB) {
    const aabbA = this.m_tree.GetFatAABB(proxyA);
    const aabbB = this.m_tree.GetFatAABB(proxyB);
    return aabbA.TestOverlap(aabbB);
  }
  GetUserData(proxy) {
    return this.m_tree.GetUserData(proxy);
  }
  GetFatAABB(proxy) {
    return this.m_tree.GetFatAABB(proxy);
  }
  GetProxyCount() {
    return this.m_proxyCount;
  }
  UpdatePairs(callback) {
    this.m_pairCount = 0;
    for (const proxy of this.m_moveBuffer) {
      this.m_queryProxy = proxy;
      const fatAABB = this.m_tree.GetFatAABB(this.m_queryProxy);
      this.m_tree.Query(this.m_queryCallback, fatAABB);
    }
    this.m_queryProxy = null;
    this.m_moveBuffer.length = 0;
    for (let i = 0; i < this.m_pairCount;) {
      const primaryPair = this.m_pairBuffer[i];
      const userDataA = this.m_tree.GetUserData(primaryPair.proxyA);
      const userDataB = this.m_tree.GetUserData(primaryPair.proxyB);
      callback(userDataA, userDataB);
      ++i;
      while (i < this.m_pairCount) {
        const pair = this.m_pairBuffer[i];
        if (pair.proxyA != primaryPair.proxyA || pair.proxyB != primaryPair.proxyB) {
          break;
        }
        ++i;
      }
    }
  }
  Query(callback, aabb) {
    this.m_tree.Query(callback, aabb);
  }
  RayCast(callback, input) {
    this.m_tree.RayCast(callback, input);
  }
  Validate() {}
  Rebalance(iterations = 0) {
    this.m_tree.Rebalance(iterations);
  }
  BufferMove(proxy) {
    this.m_moveBuffer[this.m_moveBuffer.length] = proxy;
  }
  UnBufferMove(proxy) {
    const i = this.m_moveBuffer.indexOf(proxy);
    if (i >= 0) this.m_moveBuffer.splice(i, 1);
  }
  ComparePairs(pair1, pair2) {
    return 0;
  }
};

const b2DynamicTreeNode = class {
  constructor() {
    this.aabb = new b2AABB();
    this.userData = null;
    this.parent = null;
    this.child1 = null;
    this.child2 = null;
  }
  IsLeaf() {
    return this.child1 == null;
  }
};

const b2DynamicTreePair = class {
};

const b2Point = class {
  constructor() {
    this.p = new b2Vec2();
  }
  Support(xf, vX, vY) {
    return this.p;
  }
  GetFirstVertex(xf) {
    return this.p;
  }
};

const b2RayCastInput = class {
  constructor(p1 = null, p2 = null, maxFraction = 1) {
    this.p1 = new b2Vec2();
    this.p2 = new b2Vec2();
    if (p1) this.p1.SetV(p1);
    if (p2) this.p2.SetV(p2);
    this.maxFraction = maxFraction;
  }
};

const b2RayCastOutput = class {
  constructor() {
    this.normal = new b2Vec2();
  }
};

const b2Segment = class {
  constructor() {
    this.p1 = new b2Vec2();
    this.p2 = new b2Vec2();
  }
  TestSegment(lambda, normal, segment, maxLambda = 0) {
    const s = segment.p1;
    const rX = segment.p2.x - s.x;
    const rY = segment.p2.y - s.y;
    const dX = this.p2.x - this.p1.x;
    const dY = this.p2.y - this.p1.y;
    let nX = dY;
    let nY = -dX;
    const k_slop = 100.0 * b2Settings.b2_epsilon;
    const denom = -(rX * nX + rY * nY);
    if (denom > k_slop) {
      const bX = s.x - this.p1.x;
      const bY = s.y - this.p1.y;
      let a = bX * nX + bY * nY;
      if (0.0 <= a && a <= maxLambda * denom) {
        const mu2 = (-rX * bY) + rY * bX;
        if ((-k_slop * denom) <= mu2 && mu2 <= denom * (1.0 + k_slop)) {
          a /= denom;
          const nLen = Math.sqrt(nX * nX + nY * nY);
          nX /= nLen;
          nY /= nLen;
          lambda[0] = a;
          normal.Set(nX, nY);
          return true;
        }
      }
    }
    return false;
  }
  Extend(aabb) {
    this.ExtendForward(aabb);
    this.ExtendBackward(aabb);
  }
  ExtendForward(aabb) {
    const dX = this.p2.x - this.p1.x;
    const dY = this.p2.y - this.p1.y;
    const lambda = Math.min(dX > 0 ? (aabb.upperBound.x - this.p1.x) / dX : dX < 0 ? (aabb.lowerBound.x - this.p1.x) / dX : Number.POSITIVE_INFINITY,
    dY > 0 ? (aabb.upperBound.y - this.p1.y) / dY : dY < 0 ? (aabb.lowerBound.y - this.p1.y) / dY : Number.POSITIVE_INFINITY);
    this.p2.x = this.p1.x + dX * lambda;
    this.p2.y = this.p1.y + dY * lambda;
  }
  ExtendBackward(aabb) {
    const dX = -this.p2.x + this.p1.x;
    const dY = -this.p2.y + this.p1.y;
    const lambda = Math.min(dX > 0 ? (aabb.upperBound.x - this.p2.x) / dX : dX < 0 ? (aabb.lowerBound.x - this.p2.x) / dX : Number.POSITIVE_INFINITY,
    dY > 0 ? (aabb.upperBound.y - this.p2.y) / dY : dY < 0 ? (aabb.lowerBound.y - this.p2.y) / dY : Number.POSITIVE_INFINITY);
    this.p1.x = this.p2.x + dX * lambda;
    this.p1.y = this.p2.y + dY * lambda;
  }
};

const b2SeparationFunction = class {
  static e_points = 0x01;
  static e_faceA = 0x02;
  static e_faceB = 0x04;
  static s_axis = new b2Vec2();
  constructor() {
    this.m_localPoint = new b2Vec2();
    this.m_axis = new b2Vec2();
  }
  Initialize(cache, proxyA, transformA, proxyB, transformB) {
    this.m_proxyA = proxyA;
    this.m_proxyB = proxyB;
    const count = cache.count;
    b2Settings.b2Assert(0 < count && count < 3);
    let localPointA;
    let localPointA1;
    let localPointA2;
    let localPointB;
    let localPointB1;
    let localPointB2;
    let pointAX;
    let pointAY;
    let pointBX;
    let pointBY;
    let normalX;
    let normalY;
    let tMat;
    let tVec;
    let s;
    if (count == 1) {
      this.m_type = b2SeparationFunction.e_points;
      localPointA = this.m_proxyA.GetVertex(cache.indexA[0]);
      localPointB = this.m_proxyB.GetVertex(cache.indexB[0]);
      tVec = localPointA;
      tMat = transformA.R;
      pointAX = transformA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointAY = transformA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tVec = localPointB;
      tMat = transformB.R;
      pointBX = transformB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointBY = transformB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      this.m_axis.x = pointBX - pointAX;
      this.m_axis.y = pointBY - pointAY;
      this.m_axis.Normalize();
    } else if (cache.indexB[0] == cache.indexB[1]) {
      this.m_type = b2SeparationFunction.e_faceA;
      localPointA1 = this.m_proxyA.GetVertex(cache.indexA[0]);
      localPointA2 = this.m_proxyA.GetVertex(cache.indexA[1]);
      localPointB = this.m_proxyB.GetVertex(cache.indexB[0]);
      this.m_localPoint.x = 0.5 * (localPointA1.x + localPointA2.x);
      this.m_localPoint.y = 0.5 * (localPointA1.y + localPointA2.y);
      this.m_axis.x = 1.0 * (localPointA2.y - localPointA1.y);
      this.m_axis.y = -1.0 * (localPointA2.x - localPointA1.x);
      this.m_axis.Normalize();
      tVec = this.m_axis;
      tMat = transformA.R;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = this.m_localPoint;
      tMat = transformA.R;
      pointAX = transformA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointAY = transformA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tVec = localPointB;
      tMat = transformB.R;
      pointBX = transformB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointBY = transformB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      s = (pointBX - pointAX) * normalX + (pointBY - pointAY) * normalY;
      if (s < 0.0) {
        this.m_axis.NegativeSelf();
      }
    } else if (cache.indexA[0] == cache.indexA[0]) {
      this.m_type = b2SeparationFunction.e_faceB;
      localPointB1 = this.m_proxyB.GetVertex(cache.indexB[0]);
      localPointB2 = this.m_proxyB.GetVertex(cache.indexB[1]);
      localPointA = this.m_proxyA.GetVertex(cache.indexA[0]);
      this.m_localPoint.x = 0.5 * (localPointB1.x + localPointB2.x);
      this.m_localPoint.y = 0.5 * (localPointB1.y + localPointB2.y);
      this.m_axis.x = 1.0 * (localPointB2.y - localPointB1.y);
      this.m_axis.y = -1.0 * (localPointB2.x - localPointB1.x);
      this.m_axis.Normalize();
      tVec = this.m_axis;
      tMat = transformB.R;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = this.m_localPoint;
      tMat = transformB.R;
      pointBX = transformB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointBY = transformB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tVec = localPointA;
      tMat = transformA.R;
      pointAX = transformA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      pointAY = transformA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      s = (pointAX - pointBX) * normalX + (pointAY - pointBY) * normalY;
      if (s < 0.0) {
        this.m_axis.NegativeSelf();
      }
    } else {
      localPointA1 = this.m_proxyA.GetVertex(cache.indexA[0]);
      localPointA2 = this.m_proxyA.GetVertex(cache.indexA[1]);
      localPointB1 = this.m_proxyB.GetVertex(cache.indexB[0]);
      localPointB2 = this.m_proxyB.GetVertex(cache.indexB[1]);
      const dA = b2Math.MulMV(transformA.R, b2Math.SubtractVV(localPointA2, localPointA1));
      const dB = b2Math.MulMV(transformB.R, b2Math.SubtractVV(localPointB2, localPointB1));
      const a = dA.x * dA.x + dA.y * dA.y;
      const e = dB.x * dB.x + dB.y * dB.y;
      const r = b2Math.SubtractVV(dB, dA);
      const c = dA.x * r.x + dA.y * r.y;
      const f = dB.x * r.x + dB.y * r.y;
      const b = dA.x * dB.x + dA.y * dB.y;
      const denom = a * e - b * b;
      s = 0.0;
      if (denom != 0.0) {
        s = b2Math.Clamp((b * f - c * e) / denom, 0.0, 1.0);
      }
      const t = (b * s + f) / e;
      if (t < 0.0) {
        s = b2Math.Clamp((b - c) / a, 0.0, 1.0);
      }
      localPointA = new b2Vec2();
      localPointA.x = localPointA1.x + s * (localPointA2.x - localPointA1.x);
      localPointA.y = localPointA1.y + s * (localPointA2.y - localPointA1.y);
      localPointB = new b2Vec2();
      localPointB.x = localPointB1.x + s * (localPointB2.x - localPointB1.x);
      localPointB.y = localPointB1.y + s * (localPointB2.y - localPointB1.y);
      if (s == 0.0 || s == 1.0) {
        this.m_type = b2SeparationFunction.e_faceB;
        this.m_axis = b2Math.CrossVF(b2Math.SubtractVV(localPointB2, localPointB1), 1.0);
        this.m_axis.Normalize();
        this.m_localPoint = localPointB;
        if (s < 0.0) {
          this.m_axis.NegativeSelf();
        }
      } else {
        this.m_type = b2SeparationFunction.e_faceA;
        this.m_axis = b2Math.CrossVF(b2Math.SubtractVV(localPointA2, localPointA1), 1.0);
        this.m_localPoint = localPointA;
        if (s < 0.0) {
          this.m_axis.NegativeSelf();
        }
      }
    }
  }
  Evaluate(transformA, transformB) {
    const axis = b2SeparationFunction.s_axis;
    const axisX = this.m_axis.x;
    const axisY = this.m_axis.y;
    let localPoint;
    let localPointA;
    let localPointB;
    let pointAX;
    let pointAY;
    let pointBX;
    let pointBY;
    let normalX;
    let normalY;
    let seperation;
    let tMat;
    switch (this.m_type) {
    case b2SeparationFunction.e_points:
      tMat = transformA.R;
      axis.x = axisX * tMat.col1.x + axisY * tMat.col1.y;
      axis.y = axisX * tMat.col2.x + axisY * tMat.col2.y;
      localPointA = this.m_proxyA.GetSupportVertex(axis);
      pointAX = tMat.col1.x * localPointA.x + tMat.col2.x * localPointA.y + transformA.position.x;
      pointAY = tMat.col1.y * localPointA.x + tMat.col2.y * localPointA.y + transformA.position.y;
      tMat = transformB.R;
      axis.x = -axisX * tMat.col1.x + -axisY * tMat.col1.y;
      axis.y = -axisX * tMat.col2.x + -axisY * tMat.col2.y;
      localPointB = this.m_proxyB.GetSupportVertex(axis);
      pointBX = tMat.col1.x * localPointB.x + tMat.col2.x * localPointB.y + transformB.position.x;
      pointBY = tMat.col1.y * localPointB.x + tMat.col2.y * localPointB.y + transformB.position.y;
      seperation = (pointBX - pointAX) * axisX + (pointBY - pointAY) * axisY;
      return seperation;
    case b2SeparationFunction.e_faceA:
      localPoint = this.m_localPoint;
      tMat = transformA.R;
      normalX = tMat.col1.x * axisX + tMat.col2.x * axisY;
      normalY = tMat.col1.y * axisX + tMat.col2.y * axisY;
      pointAX = tMat.col1.x * localPoint.x + tMat.col2.x * localPoint.y + transformA.position.x;
      pointAY = tMat.col1.y * localPoint.x + tMat.col2.y * localPoint.y + transformA.position.y;
      tMat = transformB.R;
      axis.x = -normalX * tMat.col1.x + -normalY * tMat.col1.y;
      axis.y = -normalX * tMat.col2.x + -normalY * tMat.col2.y;
      localPointB = this.m_proxyB.GetSupportVertex(axis);
      pointBX = tMat.col1.x * localPointB.x + tMat.col2.x * localPointB.y + transformB.position.x;
      pointBY = tMat.col1.y * localPointB.x + tMat.col2.y * localPointB.y + transformB.position.y;
      seperation = (pointBX - pointAX) * normalX + (pointBY - pointAY) * normalY;
      return seperation;
    case b2SeparationFunction.e_faceB:
      localPoint = this.m_localPoint;
      tMat = transformB.R;
      normalX = tMat.col1.x * axisX + tMat.col2.x * axisY;
      normalY = tMat.col1.y * axisX + tMat.col2.y * axisY;
      pointBX = tMat.col1.x * localPoint.x + tMat.col2.x * localPoint.y + transformB.position.x;
      pointBY = tMat.col1.y * localPoint.x + tMat.col2.y * localPoint.y + transformB.position.y;
      tMat = transformA.R;
      axis.x = -normalX * tMat.col1.x + -normalY * tMat.col1.y;
      axis.y = -normalX * tMat.col2.x + -normalY * tMat.col2.y;
      localPointA = this.m_proxyA.GetSupportVertex(axis);
      pointAX = tMat.col1.x * localPointA.x + tMat.col2.x * localPointA.y + transformA.position.x;
      pointAY = tMat.col1.y * localPointA.x + tMat.col2.y * localPointA.y + transformA.position.y;
      seperation = (pointAX - pointBX) * normalX + (pointAY - pointBY) * normalY;
      return seperation;
    default:
      b2Settings.b2Assert(false);
      return 0.0;
    }
  }
};

const b2SimplexCache = class {
  constructor() {
    this.indexA = new Array(3).fill(0);
    this.indexB = new Array(3).fill(0);
  }
};

const b2TimeOfImpact = class {
  static b2_toiCalls = 0;
  static b2_toiIters = 0;
  static b2_toiMaxIters = 0;
  static b2_toiRootIters = 0;
  static b2_toiMaxRootIters = 0;
  static s_cache = new b2SimplexCache();
  static s_distanceInput = new b2DistanceInput();
  static s_xfA = new b2Transform();
  static s_xfB = new b2Transform();
  static s_finalXfA = new b2Transform();
  static s_finalXfB = new b2Transform();
  static s_fcn = new b2SeparationFunction();
  static s_distanceOutput = new b2DistanceOutput();
  static TimeOfImpact(input) {
    ++b2TimeOfImpact.b2_toiCalls;
    const proxyA = input.proxyA;
    const proxyB = input.proxyB;
    const sweepA = input.sweepA;
    const sweepB = input.sweepB;
    b2Settings.b2Assert(sweepA.t0 == sweepB.t0);
    b2Settings.b2Assert(1.0 - sweepA.t0 > b2Settings.b2_epsilon);
    const radius = proxyA.m_radius + proxyB.m_radius;
    const tolerance = input.tolerance;
    let alpha = 0.0;
    const k_maxIterations = 1000;
    let iter = 0;
    let target = 0.0;
    const cache = b2TimeOfImpact.s_cache;
    const distanceInput = b2TimeOfImpact.s_distanceInput;
    const distanceOutput = b2TimeOfImpact.s_distanceOutput;
    const separationFunction = b2TimeOfImpact.s_fcn;
    const xfA = b2TimeOfImpact.s_xfA;
    const xfB = b2TimeOfImpact.s_xfB;
    const finalXfA = b2TimeOfImpact.s_finalXfA;
    const finalXfB = b2TimeOfImpact.s_finalXfB;
    let alphaTransformsReady = false;
    let finalTransformsReady = false;
    cache.count = 0;
    distanceInput.useRadii = false;
    distanceInput.proxyA = proxyA;
    distanceInput.proxyB = proxyB;
    distanceInput.transformA = xfA;
    distanceInput.transformB = xfB;
    for (;;) {
      if (alphaTransformsReady == false) {
        sweepA.GetTransform(xfA, alpha);
        sweepB.GetTransform(xfB, alpha);
      }
      alphaTransformsReady = false;
      b2Distance.Distance(distanceOutput, cache, distanceInput);
      if (distanceOutput.distance <= 0.0) {
        alpha = 1.0;
        break;
      }
      separationFunction.Initialize(cache, proxyA, xfA, proxyB, xfB);
      const separation = separationFunction.Evaluate(xfA, xfB);
      if (separation <= 0.0) {
        alpha = 1.0;
        break;
      }
      if (iter == 0) {
        if (separation > radius) {
          target = b2Math.Max(radius - tolerance, 0.75 * radius);
        } else {
          target = b2Math.Max(separation - tolerance, 0.02 * radius);
        }
      }
      if (separation - target < 0.5 * tolerance) {
        if (iter == 0) {
          alpha = 1.0;
          break;
        }
        break;
      }
      let newAlpha = alpha; {
        let x1 = alpha;
        let x2 = 1.0;
        let f1 = separation;
        if (finalTransformsReady == false) {
          sweepA.GetTransform(finalXfA, x2);
          sweepB.GetTransform(finalXfB, x2);
          finalTransformsReady = true;
        }
        let f2 = separationFunction.Evaluate(finalXfA, finalXfB);
        if (f2 >= target) {
          alpha = 1.0;
          break;
        }
        let rootIterCount = 0;
        for (;;) {
          let x;
          if (rootIterCount & 1) {
            x = x1 + (target - f1) * (x2 - x1) / (f2 - f1);
          } else {
            x = 0.5 * (x1 + x2);
          }
          sweepA.GetTransform(xfA, x);
          sweepB.GetTransform(xfB, x);
          const f = separationFunction.Evaluate(xfA, xfB);
          if (b2Math.Abs(f - target) < 0.025 * tolerance) {
            newAlpha = x;
            alphaTransformsReady = true;
            break;
          }
          if (f > target) {
            x1 = x;
            f1 = f;
          } else {
            x2 = x;
            f2 = f;
          }
          ++rootIterCount;
          ++b2TimeOfImpact.b2_toiRootIters;
          if (rootIterCount == 50) {
            break;
          }
        }
        b2TimeOfImpact.b2_toiMaxRootIters = b2Math.Max(b2TimeOfImpact.b2_toiMaxRootIters, rootIterCount);
      }
      if (newAlpha < (1.0 + 100.0 * b2Settings.b2_epsilon) * alpha) {
        break;
      }
      alpha = newAlpha;
      iter++;
      ++b2TimeOfImpact.b2_toiIters;
      if (iter == k_maxIterations) {
        break;
      }
    }
    b2TimeOfImpact.b2_toiMaxIters = b2Math.Max(b2TimeOfImpact.b2_toiMaxIters, iter);
    return alpha;
  }
};

const b2TOIInput = class {
  constructor() {
    this.proxyA = new b2DistanceProxy();
    this.proxyB = new b2DistanceProxy();
    this.sweepA = new b2Sweep();
    this.sweepB = new b2Sweep();
  }
};

const b2WorldManifold = class {
  constructor() {
    this.m_normal = new b2Vec2();
    this.m_points = new Array(b2Settings.b2_maxManifoldPoints);
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.m_points[i] = new b2Vec2();
    }
  }
  Initialize(manifold, xfA, radiusA = 0, xfB, radiusB = 0) {
    const pointCount = manifold.m_pointCount;
    if (pointCount == 0) {
      return;
    }
    let tVec;
    let tMat;
    let normalX;
    let normalY;
    let planePointX;
    let planePointY;
    let clipPointX;
    let clipPointY;
    let halfSeparation;
    const normal = this.m_normal;
    const points = this.m_points;
    const manifoldPoints = manifold.m_points;
    let point;
    switch (manifold.m_type) {
    case b2Manifold.e_circles: {
      tMat = xfA.R;
      tVec = manifold.m_localPoint;
      const pointAX = xfA.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      const pointAY = xfA.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tMat = xfB.R;
      tVec = manifoldPoints[0].m_localPoint;
      const pointBX = xfB.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      const pointBY = xfB.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      const dX = pointBX - pointAX;
      const dY = pointBY - pointAY;
      const d2 = dX * dX + dY * dY;
      if (d2 > b2Settings.b2_epsilon * b2Settings.b2_epsilon) {
        const d = Math.sqrt(d2);
        normalX = dX / d;
        normalY = dY / d;
      } else {
        normalX = 1;
        normalY = 0;
      }
      normal.x = normalX;
      normal.y = normalY;
      const cAX = pointAX + radiusA * normalX;
      const cAY = pointAY + radiusA * normalY;
      const cBX = pointBX - radiusB * normalX;
      const cBY = pointBY - radiusB * normalY;
      point = points[0];
      point.x = 0.5 * (cAX + cBX);
      point.y = 0.5 * (cAY + cBY);
    }
      break;
    case b2Manifold.e_faceA:
      tMat = xfA.R;
      tVec = manifold.m_localPlaneNormal;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = manifold.m_localPoint;
      planePointX = xfA.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      planePointY = xfA.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      normal.x = normalX;
      normal.y = normalY;
      tMat = xfB.R;
      for (let i = 0; i < pointCount; i++) {
        tVec = manifoldPoints[i].m_localPoint;
        clipPointX = xfB.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
        clipPointY = xfB.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
        halfSeparation = 0.5 * (radiusA - (clipPointX - planePointX) * normalX - (clipPointY - planePointY) * normalY - radiusB);
        point = points[i];
        point.x = clipPointX + halfSeparation * normalX;
        point.y = clipPointY + halfSeparation * normalY;
      }
      break;
    case b2Manifold.e_faceB:
      tMat = xfB.R;
      tVec = manifold.m_localPlaneNormal;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = manifold.m_localPoint;
      planePointX = xfB.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      planePointY = xfB.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      normal.x = -normalX;
      normal.y = -normalY;
      tMat = xfA.R;
      for (let i = 0; i < pointCount; i++) {
        tVec = manifoldPoints[i].m_localPoint;
        clipPointX = xfA.position.x + tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
        clipPointY = xfA.position.y + tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
        halfSeparation = 0.5 * (radiusB - (clipPointX - planePointX) * normalX - (clipPointY - planePointY) * normalY - radiusA);
        point = points[i];
        point.x = clipPointX + halfSeparation * normalX;
        point.y = clipPointY + halfSeparation * normalY;
      }
      break;
    }
  }
};

const b2Shape = class {
  static e_unknownShape = -1;
  static e_circleShape = 0;
  static e_polygonShape = 1;
  static e_edgeShape = 2;
  static e_shapeTypeCount = 3;
  static e_hitCollide = 1;
  static e_missCollide = 0;
  static e_startsInsideCollide = -1;
  static s_input = new b2DistanceInput();
  static {
    this.s_input.proxyA = new b2DistanceProxy();
  }
  static {
    this.s_input.proxyB = new b2DistanceProxy();
  }
  static s_simplexCache = new b2SimplexCache();
  static s_output = new b2DistanceOutput();
  constructor() {
    this.m_type = b2Shape.e_unknownShape;
    this.m_radius = b2Settings.b2_linearSlop;
  }
  Copy() {
    return null;
  }
  Set(other) {
    this.m_radius = other.m_radius;
  }
  GetType() {
    return this.m_type;
  }
  TestPoint(xf, p) {
    return false;
  }
  RayCast(output, input, transform) {
    return false;
  }
  ComputeAABB(aabb, xf) {}
  ComputeMass(massData, density) {}
  ComputeSubmergedArea(normal, offset, xf, c) {
    return 0;
  }
  static TestOverlap(shape1, transform1, shape2, transform2) {
    const input = b2Shape.s_input;
    input.proxyA.Set(shape1);
    input.proxyB.Set(shape2);
    input.transformA = transform1;
    input.transformB = transform2;
    input.useRadii = true;
    const simplexCache = b2Shape.s_simplexCache;
    simplexCache.count = 0;
    const output = b2Shape.s_output;
    b2Distance.Distance(output, simplexCache, input);
    return output.distance < 10.0 * b2Settings.b2_epsilon;
  }
};

const b2CircleShape = class extends b2Shape {
  constructor(radius = 0) {
    super();
    this.m_p = new b2Vec2();
    this.m_type = b2Shape.e_circleShape;
    this.m_radius = radius;
  }
  Copy() {
    const s = new b2CircleShape();
    s.Set(this);
    return s;
  }
  Set(other) {
    super.Set(other);
    if (other instanceof b2CircleShape) {
      const other2 = other;
      this.m_p.SetV(other2.m_p);
    }
  }
  TestPoint(transform, p) {
    const tMat = transform.R;
    let dX = transform.position.x + (tMat.col1.x * this.m_p.x + tMat.col2.x * this.m_p.y);
    let dY = transform.position.y + (tMat.col1.y * this.m_p.x + tMat.col2.y * this.m_p.y);
    dX = p.x - dX;
    dY = p.y - dY;
    return (dX * dX + dY * dY) <= this.m_radius * this.m_radius;
  }
  RayCast(output, input, transform) {
    const tMat = transform.R;
    const positionX = transform.position.x + (tMat.col1.x * this.m_p.x + tMat.col2.x * this.m_p.y);
    const positionY = transform.position.y + (tMat.col1.y * this.m_p.x + tMat.col2.y * this.m_p.y);
    const sX = input.p1.x - positionX;
    const sY = input.p1.y - positionY;
    const b = (sX * sX + sY * sY) - this.m_radius * this.m_radius;
    const rX = input.p2.x - input.p1.x;
    const rY = input.p2.y - input.p1.y;
    const c = sX * rX + sY * rY;
    const rr = rX * rX + rY * rY;
    const sigma = c * c - rr * b;
    if (sigma < 0.0 || rr < b2Settings.b2_epsilon) {
      return false;
    }
    let a = -(c + Math.sqrt(sigma));
    if (0.0 <= a && a <= input.maxFraction * rr) {
      a /= rr;
      output.fraction = a;
      output.normal.x = sX + a * rX;
      output.normal.y = sY + a * rY;
      output.normal.Normalize();
      return true;
    }
    return false;
  }
  ComputeAABB(aabb, transform) {
    const tMat = transform.R;
    const pX = transform.position.x + (tMat.col1.x * this.m_p.x + tMat.col2.x * this.m_p.y);
    const pY = transform.position.y + (tMat.col1.y * this.m_p.x + tMat.col2.y * this.m_p.y);
    aabb.lowerBound.Set(pX - this.m_radius, pY - this.m_radius);
    aabb.upperBound.Set(pX + this.m_radius, pY + this.m_radius);
  }
  ComputeMass(massData, density = 0) {
    massData.mass = density * b2Settings.b2_pi * this.m_radius * this.m_radius;
    massData.center.SetV(this.m_p);
    massData.I = massData.mass * (0.5 * this.m_radius * this.m_radius + (this.m_p.x * this.m_p.x + this.m_p.y * this.m_p.y));
  }
  ComputeSubmergedArea(normal, offset = 0, xf, c) {
    const p = b2Math.MulX(xf, this.m_p);
    const l = -(b2Math.Dot(normal, p) - offset);
    if (l < -this.m_radius + Number.MIN_VALUE) {
      return 0;
    }
    if (l > this.m_radius) {
      c.SetV(p);
      return Math.PI * this.m_radius * this.m_radius;
    }
    const r2 = this.m_radius * this.m_radius;
    const l2 = l * l;
    const area = r2 * (b2Math.Atan2(l, Math.sqrt(r2 - l2)) + Math.PI / 2) + l * Math.sqrt(r2 - l2);
    const com = -2 / 3 * ((r2 - l2) * Math.sqrt(r2 - l2)) / area;
    c.x = p.x + normal.x * com;
    c.y = p.y + normal.y * com;
    return area;
  }
  GetLocalPosition() {
    return this.m_p;
  }
  SetLocalPosition(position) {
    this.m_p.SetV(position);
  }
  GetRadius() {
    return this.m_radius;
  }
  SetRadius(radius = 0) {
    this.m_radius = radius;
  }
};

const b2EdgeChainDef = class {
  constructor() {
    this.vertexCount = 0;
    this.isALoop = true;
    this.vertices = [];
  }
};

const b2EdgeShape = class extends b2Shape {
  constructor(v1, v2) {
    super();
    this.s_supportVec = new b2Vec2();
    this.m_v1 = new b2Vec2();
    this.m_v2 = new b2Vec2();
    this.m_coreV1 = new b2Vec2();
    this.m_coreV2 = new b2Vec2();
    this.m_normal = new b2Vec2();
    this.m_direction = new b2Vec2();
    this.m_cornerDir1 = new b2Vec2();
    this.m_cornerDir2 = new b2Vec2();
    this.m_type = b2Shape.e_edgeShape;
    this.m_prevEdge = null;
    this.m_nextEdge = null;
    this.m_v1 = v1;
    this.m_v2 = v2;
    this.m_direction.Set(this.m_v2.x - this.m_v1.x, this.m_v2.y - this.m_v1.y);
    this.m_length = this.m_direction.Normalize();
    this.m_normal.Set(this.m_direction.y, -this.m_direction.x);
    this.m_coreV1.Set((-b2Settings.b2_toiSlop * (this.m_normal.x - this.m_direction.x)) + this.m_v1.x, (-b2Settings.b2_toiSlop * (this.m_normal.y - this.m_direction.y)) + this.m_v1.y);
    this.m_coreV2.Set((-b2Settings.b2_toiSlop * (this.m_normal.x + this.m_direction.x)) + this.m_v2.x, (-b2Settings.b2_toiSlop * (this.m_normal.y + this.m_direction.y)) + this.m_v2.y);
    this.m_cornerDir1 = this.m_normal;
    this.m_cornerDir2.Set(-this.m_normal.x, -this.m_normal.y);
  }
  TestPoint(transform, p) {
    return false;
  }
  RayCast(output, input, transform) {
    const rX = input.p2.x - input.p1.x;
    const rY = input.p2.y - input.p1.y;
    const tMat = transform.R;
    const v1X = transform.position.x + (tMat.col1.x * this.m_v1.x + tMat.col2.x * this.m_v1.y);
    const v1Y = transform.position.y + (tMat.col1.y * this.m_v1.x + tMat.col2.y * this.m_v1.y);
    const nX = transform.position.y + (tMat.col1.y * this.m_v2.x + tMat.col2.y * this.m_v2.y) - v1Y;
    const nY = -(transform.position.x + (tMat.col1.x * this.m_v2.x + tMat.col2.x * this.m_v2.y) - v1X);
    const k_slop = 100.0 * b2Settings.b2_epsilon;
    const denom = -(rX * nX + rY * nY);
    if (denom > k_slop) {
      const bX = input.p1.x - v1X;
      const bY = input.p1.y - v1Y;
      let a = bX * nX + bY * nY;
      if (0.0 <= a && a <= input.maxFraction * denom) {
        const mu2 = (-rX * bY) + rY * bX;
        if ((-k_slop * denom) <= mu2 && mu2 <= denom * (1.0 + k_slop)) {
          a /= denom;
          output.fraction = a;
          const nLen = Math.sqrt(nX * nX + nY * nY);
          output.normal.x = nX / nLen;
          output.normal.y = nY / nLen;
          return true;
        }
      }
    }
    return false;
  }
  ComputeAABB(aabb, transform) {
    const tMat = transform.R;
    const v1X = transform.position.x + (tMat.col1.x * this.m_v1.x + tMat.col2.x * this.m_v1.y);
    const v1Y = transform.position.y + (tMat.col1.y * this.m_v1.x + tMat.col2.y * this.m_v1.y);
    const v2X = transform.position.x + (tMat.col1.x * this.m_v2.x + tMat.col2.x * this.m_v2.y);
    const v2Y = transform.position.y + (tMat.col1.y * this.m_v2.x + tMat.col2.y * this.m_v2.y);
    if (v1X < v2X) {
      aabb.lowerBound.x = v1X;
      aabb.upperBound.x = v2X;
    } else {
      aabb.lowerBound.x = v2X;
      aabb.upperBound.x = v1X;
    }
    if (v1Y < v2Y) {
      aabb.lowerBound.y = v1Y;
      aabb.upperBound.y = v2Y;
    } else {
      aabb.lowerBound.y = v2Y;
      aabb.upperBound.y = v1Y;
    }
  }
  ComputeMass(massData, density) {
    massData.mass = 0;
    massData.center.SetV(this.m_v1);
    massData.I = 0;
  }
  ComputeSubmergedArea(normal, offset = 0, xf, c) {
    const v0 = new b2Vec2(normal.x * offset, normal.y * offset);
    const v1 = b2Math.MulX(xf, this.m_v1);
    const v2 = b2Math.MulX(xf, this.m_v2);
    const d1 = b2Math.Dot(normal, v1) - offset;
    const d2 = b2Math.Dot(normal, v2) - offset;
    if (d1 > 0) {
      if (d2 > 0) {
        return 0;
      } else {
        v1.x = (-d2 / (d1 - d2) * v1.x) + d1 / (d1 - d2) * v2.x;
        v1.y = (-d2 / (d1 - d2) * v1.y) + d1 / (d1 - d2) * v2.y;
      }
    } else {
      if (d2 > 0) {
        v2.x = (-d2 / (d1 - d2) * v1.x) + d1 / (d1 - d2) * v2.x;
        v2.y = (-d2 / (d1 - d2) * v1.y) + d1 / (d1 - d2) * v2.y;
      }
    }
    c.x = (v0.x + v1.x + v2.x) / 3;
    c.y = (v0.y + v1.y + v2.y) / 3;
    return 0.5 * ((v1.x - v0.x) * (v2.y - v0.y) - (v1.y - v0.y) * (v2.x - v0.x));
  }
  GetLength() {
    return this.m_length;
  }
  GetVertex1() {
    return this.m_v1;
  }
  GetVertex2() {
    return this.m_v2;
  }
  GetCoreVertex1() {
    return this.m_coreV1;
  }
  GetCoreVertex2() {
    return this.m_coreV2;
  }
  GetNormalVector() {
    return this.m_normal;
  }
  GetDirectionVector() {
    return this.m_direction;
  }
  GetCorner1Vector() {
    return this.m_cornerDir1;
  }
  GetCorner2Vector() {
    return this.m_cornerDir2;
  }
  Corner1IsConvex() {
    return this.m_cornerConvex1;
  }
  Corner2IsConvex() {
    return this.m_cornerConvex2;
  }
  GetFirstVertex(xf) {
    const tMat = xf.R;
    return new b2Vec2(xf.position.x + (tMat.col1.x * this.m_coreV1.x + tMat.col2.x * this.m_coreV1.y), xf.position.y + (tMat.col1.y * this.m_coreV1.x + tMat.col2.y * this.m_coreV1.y));
  }
  GetNextEdge() {
    return this.m_nextEdge;
  }
  GetPrevEdge() {
    return this.m_prevEdge;
  }
  Support(xf, dX = 0, dY = 0) {
    const tMat = xf.R;
    const v1X = xf.position.x + (tMat.col1.x * this.m_coreV1.x + tMat.col2.x * this.m_coreV1.y);
    const v1Y = xf.position.y + (tMat.col1.y * this.m_coreV1.x + tMat.col2.y * this.m_coreV1.y);
    const v2X = xf.position.x + (tMat.col1.x * this.m_coreV2.x + tMat.col2.x * this.m_coreV2.y);
    const v2Y = xf.position.y + (tMat.col1.y * this.m_coreV2.x + tMat.col2.y * this.m_coreV2.y);
    if ((v1X * dX + v1Y * dY) > (v2X * dX + v2Y * dY)) {
      this.s_supportVec.x = v1X;
      this.s_supportVec.y = v1Y;
    } else {
      this.s_supportVec.x = v2X;
      this.s_supportVec.y = v2Y;
    }
    return this.s_supportVec;
  }
  SetPrevEdge(edge, core, cornerDir, convex) {
    this.m_prevEdge = edge;
    this.m_coreV1 = core;
    this.m_cornerDir1 = cornerDir;
    this.m_cornerConvex1 = convex;
  }
  SetNextEdge(edge, core, cornerDir, convex) {
    this.m_nextEdge = edge;
    this.m_coreV2 = core;
    this.m_cornerDir2 = cornerDir;
    this.m_cornerConvex2 = convex;
  }
};

const b2MassData = class {
  constructor() {
    this.mass = 0.0;
    this.center = new b2Vec2(0, 0);
    this.I = 0.0;
  }
};

const b2PolygonShape = class extends b2Shape {
  static s_mat = new b2Mat22();
  constructor() {
    super();
    this.m_type = b2Shape.e_polygonShape;
    this.m_radius = b2Settings.b2_polygonRadius;
    this.m_centroid = new b2Vec2();
    this.m_vertices = [];
    this.m_normals = [];
  }
  Copy() {
    const s = new b2PolygonShape();
    s.Set(this);
    return s;
  }
  Set(other) {
    super.Set(other);
    if (other instanceof b2PolygonShape) {
      const other2 = other;
      this.m_centroid.SetV(other2.m_centroid);
      this.m_vertexCount = other2.m_vertexCount;
      this.Reserve(this.m_vertexCount);
      for (let i = 0; i < this.m_vertexCount; i++) {
        this.m_vertices[i].SetV(other2.m_vertices[i]);
        this.m_normals[i].SetV(other2.m_normals[i]);
      }
    }
  }
  SetAsArray(vertices, vertexCount = 0) {
    const v = [];
    let i,
      tVec;
    for (i = 0; i < vertices.length; ++i) {
      tVec = vertices[i];
      v.push(tVec);
    }
    this.SetAsVector(v, vertexCount);
  }
  static AsArray(vertices, vertexCount = 0) {
    const polygonShape = new b2PolygonShape();
    polygonShape.SetAsArray(vertices, vertexCount);
    return polygonShape;
  }
  SetAsVector(vertices, vertexCount = 0) {
    if (vertexCount == 0) vertexCount = vertices.length;
    b2Settings.b2Assert(2 <= vertexCount);
    this.m_vertexCount = vertexCount;
    this.Reserve(vertexCount);
    for (let i = 0; i < this.m_vertexCount; i++) {
      this.m_vertices[i].SetV(vertices[i]);
    }
    for (let i = 0; i < this.m_vertexCount; ++i) {
      const i1 = i;
      const i2 = i + 1 < this.m_vertexCount ? i + 1 : 0;
      const edge = b2Math.SubtractVV(this.m_vertices[i2], this.m_vertices[i1]);
      b2Settings.b2Assert(edge.LengthSquared() > b2Settings.b2_epsilon * b2Settings.b2_epsilon);
      this.m_normals[i].SetV(b2Math.CrossVF(edge, 1.0));
      this.m_normals[i].Normalize();
    }
    this.m_centroid = b2PolygonShape.ComputeCentroid(this.m_vertices, this.m_vertexCount);
  }
  static AsVector(vertices, vertexCount = 0) {
    const polygonShape = new b2PolygonShape();
    polygonShape.SetAsVector(vertices, vertexCount);
    return polygonShape;
  }
  SetAsBox(hx = 0, hy = 0) {
    this.m_vertexCount = 4;
    this.Reserve(4);
    this.m_vertices[0].Set(-hx, -hy);
    this.m_vertices[1].Set(hx, -hy);
    this.m_vertices[2].Set(hx, hy);
    this.m_vertices[3].Set(-hx, hy);
    this.m_normals[0].Set(0.0, -1.0);
    this.m_normals[1].Set(1.0, 0.0);
    this.m_normals[2].Set(0.0, 1.0);
    this.m_normals[3].Set(-1.0, 0.0);
    this.m_centroid.SetZero();
  }
  static AsBox(hx = 0, hy = 0) {
    const polygonShape = new b2PolygonShape();
    polygonShape.SetAsBox(hx, hy);
    return polygonShape;
  }
  SetAsOrientedBox(hx = 0, hy = 0, center = null, angle = 0.0) {
    this.m_vertexCount = 4;
    this.Reserve(4);
    this.m_vertices[0].Set(-hx, -hy);
    this.m_vertices[1].Set(hx, -hy);
    this.m_vertices[2].Set(hx, hy);
    this.m_vertices[3].Set(-hx, hy);
    this.m_normals[0].Set(0.0, -1.0);
    this.m_normals[1].Set(1.0, 0.0);
    this.m_normals[2].Set(0.0, 1.0);
    this.m_normals[3].Set(-1.0, 0.0);
    this.m_centroid.SetV(center);
    const xf = new b2Transform();
    xf.position = center;
    xf.R.Set(angle);
    for (let i = 0; i < this.m_vertexCount; ++i) {
      this.m_vertices[i] = b2Math.MulX(xf, this.m_vertices[i]);
      this.m_normals[i] = b2Math.MulMV(xf.R, this.m_normals[i]);
    }
  }
  static AsOrientedBox(hx = 0, hy = 0, center = null, angle = 0.0) {
    const polygonShape = new b2PolygonShape();
    polygonShape.SetAsOrientedBox(hx, hy, center, angle);
    return polygonShape;
  }
  SetAsEdge(v1, v2) {
    this.m_vertexCount = 2;
    this.Reserve(2);
    this.m_vertices[0].SetV(v1);
    this.m_vertices[1].SetV(v2);
    this.m_centroid.x = 0.5 * (v1.x + v2.x);
    this.m_centroid.y = 0.5 * (v1.y + v2.y);
    this.m_normals[0] = b2Math.CrossVF(b2Math.SubtractVV(v2, v1), 1.0);
    this.m_normals[0].Normalize();
    this.m_normals[1].x = -this.m_normals[0].x;
    this.m_normals[1].y = -this.m_normals[0].y;
  }
  static AsEdge(v1, v2) {
    const polygonShape = new b2PolygonShape();
    polygonShape.SetAsEdge(v1, v2);
    return polygonShape;
  }
  TestPoint(xf, p) {
    let tVec;
    const tMat = xf.R;
    let tX = p.x - xf.position.x;
    let tY = p.y - xf.position.y;
    const pLocalX = tX * tMat.col1.x + tY * tMat.col1.y;
    const pLocalY = tX * tMat.col2.x + tY * tMat.col2.y;
    for (let i = 0; i < this.m_vertexCount; ++i) {
      tVec = this.m_vertices[i];
      tX = pLocalX - tVec.x;
      tY = pLocalY - tVec.y;
      tVec = this.m_normals[i];
      const dot = tVec.x * tX + tVec.y * tY;
      if (dot > 0.0) {
        return false;
      }
    }
    return true;
  }
  RayCast(output, input, transform) {
    let lower = 0.0;
    let upper = input.maxFraction;
    let tVec;
    let tX = input.p1.x - transform.position.x;
    let tY = input.p1.y - transform.position.y;
    let tMat = transform.R;
    const p1X = tX * tMat.col1.x + tY * tMat.col1.y;
    const p1Y = tX * tMat.col2.x + tY * tMat.col2.y;
    tX = input.p2.x - transform.position.x;
    tY = input.p2.y - transform.position.y;
    tMat = transform.R;
    const p2X = tX * tMat.col1.x + tY * tMat.col1.y;
    const p2Y = tX * tMat.col2.x + tY * tMat.col2.y;
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    let index = -1;
    for (let i = 0; i < this.m_vertexCount; ++i) {
      tVec = this.m_vertices[i];
      tX = tVec.x - p1X;
      tY = tVec.y - p1Y;
      tVec = this.m_normals[i];
      const numerator = tVec.x * tX + tVec.y * tY;
      const denominator = tVec.x * dX + tVec.y * dY;
      if (denominator == 0.0) {
        if (numerator < 0.0) {
          return false;
        }
      } else {
        if (denominator < 0.0 && numerator < lower * denominator) {
          lower = numerator / denominator;
          index = i;
        } else if (denominator > 0.0 && numerator < upper * denominator) {
          upper = numerator / denominator;
        }
      }
      if (upper < lower) {
        return false;
      }
    }
    if (index >= 0) {
      output.fraction = lower;
      tMat = transform.R;
      tVec = this.m_normals[index];
      output.normal.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      output.normal.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      return true;
    }
    return false;
  }
  ComputeAABB(aabb, xf) {
    const tMat = xf.R;
    let tVec = this.m_vertices[0];
    let lowerX = xf.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    let lowerY = xf.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    let upperX = lowerX;
    let upperY = lowerY;
    for (let i = 1; i < this.m_vertexCount; ++i) {
      tVec = this.m_vertices[i];
      const vX = xf.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      const vY = xf.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      lowerX = lowerX < vX ? lowerX : vX;
      lowerY = lowerY < vY ? lowerY : vY;
      upperX = upperX > vX ? upperX : vX;
      upperY = upperY > vY ? upperY : vY;
    }
    aabb.lowerBound.x = lowerX - this.m_radius;
    aabb.lowerBound.y = lowerY - this.m_radius;
    aabb.upperBound.x = upperX + this.m_radius;
    aabb.upperBound.y = upperY + this.m_radius;
  }
  ComputeMass(massData, density = 0) {
    if (this.m_vertexCount == 2) {
      massData.center.x = 0.5 * (this.m_vertices[0].x + this.m_vertices[1].x);
      massData.center.y = 0.5 * (this.m_vertices[0].y + this.m_vertices[1].y);
      massData.mass = 0.0;
      massData.I = 0.0;
      return;
    }
    let centerX = 0.0;
    let centerY = 0.0;
    let area = 0.0;
    let I = 0.0;
    const p1X = 0.0;
    const p1Y = 0.0;
    const k_inv3 = 1.0 / 3.0;
    for (let i = 0; i < this.m_vertexCount; ++i) {
      const p2 = this.m_vertices[i];
      const p3 = i + 1 < this.m_vertexCount ? this.m_vertices[i + 1] : this.m_vertices[0];
      const e1X = p2.x - p1X;
      const e1Y = p2.y - p1Y;
      const e2X = p3.x - p1X;
      const e2Y = p3.y - p1Y;
      const D = e1X * e2Y - e1Y * e2X;
      const triangleArea = 0.5 * D;area += triangleArea;
      centerX += triangleArea * k_inv3 * (p1X + p2.x + p3.x);
      centerY += triangleArea * k_inv3 * (p1Y + p2.y + p3.y);
      const px = p1X;
      const py = p1Y;
      const ex1 = e1X;
      const ey1 = e1Y;
      const ex2 = e2X;
      const ey2 = e2Y;
      const intx2 = k_inv3 * (0.25 * (ex1 * ex1 + ex2 * ex1 + ex2 * ex2) + (px * ex1 + px * ex2)) + 0.5 * px * px;
      const inty2 = k_inv3 * (0.25 * (ey1 * ey1 + ey2 * ey1 + ey2 * ey2) + (py * ey1 + py * ey2)) + 0.5 * py * py;I += D * (intx2 + inty2);
    }
    massData.mass = density * area;
    centerX *= 1.0 / area;
    centerY *= 1.0 / area;
    massData.center.Set(centerX, centerY);
    massData.I = density * I;
  }
  ComputeSubmergedArea(normal, offset = 0, xf, c) {
    const normalL = b2Math.MulTMV(xf.R, normal);
    const offsetL = offset - b2Math.Dot(normal, xf.position);
    const depths = [];
    let diveCount = 0;
    let intoIndex = -1;
    let outoIndex = -1;
    let lastSubmerged = false;
    let i;
    for (i = 0; i < this.m_vertexCount; ++i) {
      depths[i] = b2Math.Dot(normalL, this.m_vertices[i]) - offsetL;
      const isSubmerged = depths[i] < -Number.MIN_VALUE;
      if (i > 0) {
        if (isSubmerged) {
          if (!lastSubmerged) {
            intoIndex = i - 1;
            diveCount++;
          }
        } else {
          if (lastSubmerged) {
            outoIndex = i - 1;
            diveCount++;
          }
        }
      }
      lastSubmerged = isSubmerged;
    }
    switch (diveCount) {
    case 0:
      if (lastSubmerged) {
        const md = new b2MassData();
        this.ComputeMass(md, 1);
        c.SetV(b2Math.MulX(xf, md.center));
        return md.mass;
      } else {
        return 0;
      }
    case 1:
      if (intoIndex == -1) {
        intoIndex = this.m_vertexCount - 1;
      } else {
        outoIndex = this.m_vertexCount - 1;
      }
      break;
    }
    const intoIndex2 = (intoIndex + 1) % this.m_vertexCount;
    const outoIndex2 = (outoIndex + 1) % this.m_vertexCount;
    const intoLamdda = (0 - depths[intoIndex]) / (depths[intoIndex2] - depths[intoIndex]);
    const outoLamdda = (0 - depths[outoIndex]) / (depths[outoIndex2] - depths[outoIndex]);
    const intoVec = new b2Vec2(this.m_vertices[intoIndex].x * (1 - intoLamdda) + this.m_vertices[intoIndex2].x * intoLamdda, this.m_vertices[intoIndex].y * (1 - intoLamdda) + this.m_vertices[intoIndex2].y * intoLamdda);
    const outoVec = new b2Vec2(this.m_vertices[outoIndex].x * (1 - outoLamdda) + this.m_vertices[outoIndex2].x * outoLamdda, this.m_vertices[outoIndex].y * (1 - outoLamdda) + this.m_vertices[outoIndex2].y * outoLamdda);
    let area = 0;
    const center = new b2Vec2();
    let p2 = this.m_vertices[intoIndex2];
    let p3;
    i = intoIndex2;
    while (i != outoIndex2) {
      i = (i + 1) % this.m_vertexCount;
      if (i == outoIndex2) p3 = outoVec;
      else p3 = this.m_vertices[i];
      const triangleArea = 0.5 * ((p2.x - intoVec.x) * (p3.y - intoVec.y) - (p2.y - intoVec.y) * (p3.x - intoVec.x));
      area += triangleArea;
      center.x += triangleArea * (intoVec.x + p2.x + p3.x) / 3;
      center.y += triangleArea * (intoVec.y + p2.y + p3.y) / 3;
      p2 = p3;
    }
    center.Multiply(1 / area);
    c.SetV(b2Math.MulX(xf, center));
    return area;
  }
  GetVertexCount() {
    return this.m_vertexCount;
  }
  GetVertices() {
    return this.m_vertices;
  }
  GetNormals() {
    return this.m_normals;
  }
  GetSupport(d) {
    let bestIndex = 0;
    let bestValue = this.m_vertices[0].x * d.x + this.m_vertices[0].y * d.y;
    for (let i = 1; i < this.m_vertexCount; ++i) {
      const value = this.m_vertices[i].x * d.x + this.m_vertices[i].y * d.y;
      if (value > bestValue) {
        bestIndex = i;
        bestValue = value;
      }
    }
    return bestIndex;
  }
  GetSupportVertex(d) {
    let bestIndex = 0;
    let bestValue = this.m_vertices[0].x * d.x + this.m_vertices[0].y * d.y;
    for (let i = 1; i < this.m_vertexCount; ++i) {
      const value = this.m_vertices[i].x * d.x + this.m_vertices[i].y * d.y;
      if (value > bestValue) {
        bestIndex = i;
        bestValue = value;
      }
    }
    return this.m_vertices[bestIndex];
  }
  Validate() {
    return false;
  }
  Reserve(count = 0) {
    for (let i = this.m_vertices.length; i < count; i++) {
      this.m_vertices[i] = new b2Vec2();
      this.m_normals[i] = new b2Vec2();
    }
  }
  static ComputeCentroid(vs, count = 0) {
    const c = new b2Vec2();
    let area = 0.0;
    const p1X = 0.0;
    const p1Y = 0.0;
    const inv3 = 1.0 / 3.0;
    for (let i = 0; i < count; ++i) {
      const p2 = vs[i];
      const p3 = i + 1 < count ? vs[i + 1] : vs[0];
      const e1X = p2.x - p1X;
      const e1Y = p2.y - p1Y;
      const e2X = p3.x - p1X;
      const e2Y = p3.y - p1Y;
      const D = e1X * e2Y - e1Y * e2X;
      const triangleArea = 0.5 * D;area += triangleArea;
      c.x += triangleArea * inv3 * (p1X + p2.x + p3.x);
      c.y += triangleArea * inv3 * (p1Y + p2.y + p3.y);
    }
    c.x *= 1.0 / area;
    c.y *= 1.0 / area;
    return c;
  }
  static ComputeOBB(obb, vs, count = 0) {
    const p = new Array(count + 1);
    for (let i = 0; i < count; ++i) {
      p[i] = vs[i];
    }
    p[count] = p[0];
    let minArea = Number.MAX_VALUE;
    for (let i = 1; i <= count; ++i) {
      const root = p[i - 1];
      let uxX = p[i].x - root.x;
      let uxY = p[i].y - root.y;
      const length = Math.sqrt(uxX * uxX + uxY * uxY);
      uxX /= length;
      uxY /= length;
      const uyX = -uxY;
      const uyY = uxX;
      let lowerX = Number.MAX_VALUE;
      let lowerY = Number.MAX_VALUE;
      let upperX = -Number.MAX_VALUE;
      let upperY = -Number.MAX_VALUE;
      for (let j = 0; j < count; ++j) {
        const dX = p[j].x - root.x;
        const dY = p[j].y - root.y;
        const rX = uxX * dX + uxY * dY;
        const rY = uyX * dX + uyY * dY;
        if (rX < lowerX) lowerX = rX;
        if (rY < lowerY) lowerY = rY;
        if (rX > upperX) upperX = rX;
        if (rY > upperY) upperY = rY;
      }
      const area = (upperX - lowerX) * (upperY - lowerY);
      if (area < 0.95 * minArea) {
        minArea = area;
        obb.R.col1.x = uxX;
        obb.R.col1.y = uxY;
        obb.R.col2.x = uyX;
        obb.R.col2.y = uyY;
        const centerX = 0.5 * (lowerX + upperX);
        const centerY = 0.5 * (lowerY + upperY);
        const tMat = obb.R;
        obb.center.x = root.x + (tMat.col1.x * centerX + tMat.col2.x * centerY);
        obb.center.y = root.y + (tMat.col1.y * centerX + tMat.col2.y * centerY);
        obb.extents.x = 0.5 * (upperX - lowerX);
        obb.extents.y = 0.5 * (upperY - lowerY);
      }
    }
  }
};

const b2Body = class {
  static e_islandFlag = 0x0001;
  static e_awakeFlag = 0x0002;
  static e_allowSleepFlag = 0x0004;
  static e_bulletFlag = 0x0008;
  static e_fixedRotationFlag = 0x0010;
  static e_activeFlag = 0x0020;
  static b2_staticBody = 0;
  static b2_kinematicBody = 1;
  static b2_dynamicBody = 2;
  static s_xf1 = new b2Transform();
  constructor(bd, world) {
    this.m_xf = new b2Transform();
    this.m_sweep = new b2Sweep();
    this.m_linearVelocity = new b2Vec2();
    this.m_force = new b2Vec2();
    this.m_islandIndex = 0;
    this.m_flags = 0;
    if (bd.bullet) {
      this.m_flags |= b2Body.e_bulletFlag;
    }
    if (bd.fixedRotation) {
      this.m_flags |= b2Body.e_fixedRotationFlag;
    }
    if (bd.allowSleep) {
      this.m_flags |= b2Body.e_allowSleepFlag;
    }
    if (bd.awake) {
      this.m_flags |= b2Body.e_awakeFlag;
    }
    if (bd.active) {
      this.m_flags |= b2Body.e_activeFlag;
    }
    this.m_world = world;
    this.m_xf.position.SetV(bd.position);
    this.m_xf.R.Set(bd.angle);
    this.m_sweep.localCenter.SetZero();
    this.m_sweep.t0 = 1.0;
    this.m_sweep.a0 = this.m_sweep.a = bd.angle;
    const tMat = this.m_xf.R;
    const tVec = this.m_sweep.localCenter;
    this.m_sweep.c.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_sweep.c.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_sweep.c.x += this.m_xf.position.x;
    this.m_sweep.c.y += this.m_xf.position.y;
    this.m_sweep.c0.SetV(this.m_sweep.c);
    this.m_jointList = null;
    this.m_controllerList = null;
    this.m_contactList = null;
    this.m_controllerCount = 0;
    this.m_prev = null;
    this.m_next = null;
    this.m_linearVelocity.SetV(bd.linearVelocity);
    this.m_angularVelocity = bd.angularVelocity;
    this.m_linearDamping = bd.linearDamping;
    this.m_angularDamping = bd.angularDamping;
    this.m_force.Set(0.0, 0.0);
    this.m_torque = 0.0;
    this.m_sleepTime = 0.0;
    this.m_type = bd.type;
    if (this.m_type == b2Body.b2_dynamicBody) {
      this.m_mass = 1.0;
      this.m_invMass = 1.0;
    } else {
      this.m_mass = 0.0;
      this.m_invMass = 0.0;
    }
    this.m_I = 0.0;
    this.m_invI = 0.0;
    this.m_inertiaScale = bd.inertiaScale;
    this.m_userData = bd.userData;
    this.m_fixtureList = null;
    this.m_fixtureCount = 0;
  }
  connectEdges(s1, s2, angle1 = 0) {
    const angle2 = b2Math.Atan2(s2.GetDirectionVector().y, s2.GetDirectionVector().x);
    const sinCos = b2Math.s_sinCos;
    b2Math.SinCos((angle2 - angle1) * 0.5, sinCos);
    const coreOffset = sinCos[0] / sinCos[1];
    let core = b2Math.MulFV(coreOffset, s2.GetDirectionVector());
    core = b2Math.SubtractVV(core, s2.GetNormalVector());
    core = b2Math.MulFV(b2Settings.b2_toiSlop, core);
    core = b2Math.AddVV(core, s2.GetVertex1());
    const cornerDir = b2Math.AddVV(s1.GetDirectionVector(), s2.GetDirectionVector());
    cornerDir.Normalize();
    const convex = b2Math.Dot(s1.GetDirectionVector(), s2.GetNormalVector()) > 0.0;
    s1.SetNextEdge(s2, core, cornerDir, convex);
    s2.SetPrevEdge(s1, core, cornerDir, convex);
    return angle2;
  }
  CreateFixture(def) {
    if (this.m_world.IsLocked() == true) {
      return null;
    }
    const fixture = new b2Fixture();
    fixture.Create(this, this.m_xf, def);
    if (this.m_flags & b2Body.e_activeFlag) {
      const broadPhase = this.m_world.m_contactManager.m_broadPhase;
      fixture.CreateProxy(broadPhase, this.m_xf);
    }
    fixture.m_next = this.m_fixtureList;
    this.m_fixtureList = fixture;
    ++this.m_fixtureCount;
    fixture.m_body = this;
    if (fixture.m_density > 0.0) {
      this.ResetMassData();
    }
    this.m_world.m_flags |= b2World.e_newFixture;
    return fixture;
  }
  CreateFixture2(shape, density = 0.0) {
    const def = new b2FixtureDef();
    def.shape = shape;
    def.density = density;
    return this.CreateFixture(def);
  }
  DestroyFixture(fixture) {
    if (this.m_world.IsLocked() == true) {
      return;
    }
    let node = this.m_fixtureList;
    let ppF = null;
    while (node != null) {
      if (node == fixture) {
        if (ppF) ppF.m_next = fixture.m_next;
        else this.m_fixtureList = fixture.m_next;
        break;
      }
      ppF = node;
      node = node.m_next;
    }
    let edge = this.m_contactList;
    while (edge) {
      const c = edge.contact;
      edge = edge.next;
      const fixtureA = c.GetFixtureA();
      const fixtureB = c.GetFixtureB();
      if (fixture == fixtureA || fixture == fixtureB) {
        this.m_world.m_contactManager.Destroy(c);
      }
    }
    if (this.m_flags & b2Body.e_activeFlag) {
      const broadPhase = this.m_world.m_contactManager.m_broadPhase;
      fixture.DestroyProxy(broadPhase);
    }
    fixture.Destroy();
    fixture.m_body = null;
    fixture.m_next = null;
    --this.m_fixtureCount;
    this.ResetMassData();
  }
  SetPositionAndAngle(position, angle = 0) {
    if (this.m_world.IsLocked() == true) {
      return;
    }
    this.m_xf.R.Set(angle);
    this.m_xf.position.SetV(position);
    const tMat = this.m_xf.R;
    const tVec = this.m_sweep.localCenter;
    this.m_sweep.c.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_sweep.c.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_sweep.c.x += this.m_xf.position.x;
    this.m_sweep.c.y += this.m_xf.position.y;
    this.m_sweep.c0.SetV(this.m_sweep.c);
    this.m_sweep.a0 = this.m_sweep.a = angle;
    const broadPhase = this.m_world.m_contactManager.m_broadPhase;
    for (let f = this.m_fixtureList; f; f = f.m_next) {
      f.Synchronize(broadPhase, this.m_xf, this.m_xf);
    }
    this.m_world.m_contactManager.FindNewContacts();
  }
  SetTransform(xf) {
    this.SetPositionAndAngle(xf.position, xf.GetAngle());
  }
  GetTransform() {
    return this.m_xf;
  }
  GetPosition() {
    return this.m_xf.position;
  }
  SetPosition(position) {
    this.SetPositionAndAngle(position, this.GetAngle());
  }
  GetAngle() {
    return this.m_sweep.a;
  }
  SetAngle(angle = 0) {
    this.SetPositionAndAngle(this.GetPosition(), angle);
  }
  GetWorldCenter() {
    return this.m_sweep.c;
  }
  GetLocalCenter() {
    return this.m_sweep.localCenter;
  }
  SetLinearVelocity(v) {
    if (this.m_type == b2Body.b2_staticBody) {
      return;
    }
    this.m_linearVelocity.SetV(v);
  }
  GetLinearVelocity() {
    return this.m_linearVelocity;
  }
  SetAngularVelocity(omega = 0) {
    if (this.m_type == b2Body.b2_staticBody) {
      return;
    }
    this.m_angularVelocity = omega;
  }
  GetAngularVelocity() {
    return this.m_angularVelocity;
  }
  GetDefinition() {
    const bd = new b2BodyDef();
    bd.type = this.GetType();
    bd.allowSleep = (this.m_flags & b2Body.e_allowSleepFlag) == b2Body.e_allowSleepFlag;
    bd.angle = this.GetAngle();
    bd.angularDamping = this.m_angularDamping;
    bd.angularVelocity = this.m_angularVelocity;
    bd.fixedRotation = (this.m_flags & b2Body.e_fixedRotationFlag) == b2Body.e_fixedRotationFlag;
    bd.bullet = (this.m_flags & b2Body.e_bulletFlag) == b2Body.e_bulletFlag;
    bd.awake = (this.m_flags & b2Body.e_awakeFlag) == b2Body.e_awakeFlag;
    bd.linearDamping = this.m_linearDamping;
    bd.linearVelocity.SetV(this.GetLinearVelocity());
    bd.position = this.GetPosition();
    bd.userData = this.GetUserData();
    return bd;
  }
  ApplyForce(force, point) {
    if (this.m_type != b2Body.b2_dynamicBody) {
      return;
    }
    if (this.IsAwake() == false) {
      this.SetAwake(true);
    }
    this.m_force.x += force.x;
    this.m_force.y += force.y;
    this.m_torque += ((point.x - this.m_sweep.c.x) * force.y - (point.y - this.m_sweep.c.y) * force.x);
  }
  ApplyTorque(torque = 0) {
    if (this.m_type != b2Body.b2_dynamicBody) {
      return;
    }
    if (this.IsAwake() == false) {
      this.SetAwake(true);
    }
    this.m_torque += torque;
  }
  ApplyImpulse(impulse, point) {
    if (this.m_type != b2Body.b2_dynamicBody) {
      return;
    }
    if (this.IsAwake() == false) {
      this.SetAwake(true);
    }
    this.m_linearVelocity.x += this.m_invMass * impulse.x;
    this.m_linearVelocity.y += this.m_invMass * impulse.y;
    this.m_angularVelocity += this.m_invI * ((point.x - this.m_sweep.c.x) * impulse.y - (point.y - this.m_sweep.c.y) * impulse.x);
  }
  Split(callback) {
    const linearVelocity = this.GetLinearVelocity().Copy();
    const angularVelocity = this.GetAngularVelocity();
    const center = this.GetWorldCenter();
    const body2 = this.m_world.CreateBody(this.GetDefinition());
    let prev;
    for (let f = this.m_fixtureList; f;) {
      if (callback(f)) {
        const next = f.m_next;
        if (prev) {
          prev.m_next = next;
        } else {
          this.m_fixtureList = next;
        }
        this.m_fixtureCount--;
        f.m_next = body2.m_fixtureList;
        body2.m_fixtureList = f;
        body2.m_fixtureCount++;
        f.m_body = body2;
        f = next;
      } else {
        prev = f;
        f = f.m_next;
      }
    }
    this.ResetMassData();
    body2.ResetMassData();
    const center1 = this.GetWorldCenter();
    const center2 = body2.GetWorldCenter();
    const velocity1 = b2Math.AddVV(linearVelocity, b2Math.CrossFV(angularVelocity, b2Math.SubtractVV(center1, center)));
    const velocity2 = b2Math.AddVV(linearVelocity, b2Math.CrossFV(angularVelocity, b2Math.SubtractVV(center2, center)));
    this.SetLinearVelocity(velocity1);
    body2.SetLinearVelocity(velocity2);
    this.SetAngularVelocity(angularVelocity);
    body2.SetAngularVelocity(angularVelocity);
    this.SynchronizeFixtures();
    body2.SynchronizeFixtures();
    return body2;
  }
  Merge(other) {
    for (let f = other.m_fixtureList; f;) {
      const next = f.m_next;
      other.m_fixtureCount--;
      f.m_next = this.m_fixtureList;
      this.m_fixtureList = f;
      this.m_fixtureCount++;
      f.m_body = this;
      f = next;
    }
    other.m_fixtureList = null;
    this.ResetMassData();
    this.SynchronizeFixtures();
  }
  GetMass() {
    return this.m_mass;
  }
  GetInertia() {
    return this.m_I;
  }
  GetMassData(data) {
    data.mass = this.m_mass;
    data.I = this.m_I;
    data.center.SetV(this.m_sweep.localCenter);
  }
  SetMassData(massData) {
    b2Settings.b2Assert(this.m_world.IsLocked() == false);
    if (this.m_world.IsLocked() == true) {
      return;
    }
    if (this.m_type != b2Body.b2_dynamicBody) {
      return;
    }
    this.m_invMass = 0.0;
    this.m_I = 0.0;
    this.m_invI = 0.0;
    this.m_mass = massData.mass;
    if (this.m_mass <= 0.0) {
      this.m_mass = 1.0;
    }
    this.m_invMass = 1.0 / this.m_mass;
    if (massData.I > 0.0 && (this.m_flags & b2Body.e_fixedRotationFlag) == 0) {
      this.m_I = massData.I - this.m_mass * (massData.center.x * massData.center.x + massData.center.y * massData.center.y);
      this.m_invI = 1.0 / this.m_I;
    }
    const oldCenter = this.m_sweep.c.Copy();
    this.m_sweep.localCenter.SetV(massData.center);
    this.m_sweep.c0.SetV(b2Math.MulX(this.m_xf, this.m_sweep.localCenter));
    this.m_sweep.c.SetV(this.m_sweep.c0);
    this.m_linearVelocity.x += this.m_angularVelocity * -(this.m_sweep.c.y - oldCenter.y);
    this.m_linearVelocity.y += this.m_angularVelocity * (+(this.m_sweep.c.x - oldCenter.x));
  }
  ResetMassData() {
    this.m_mass = 0.0;
    this.m_invMass = 0.0;
    this.m_I = 0.0;
    this.m_invI = 0.0;
    this.m_sweep.localCenter.SetZero();
    if (this.m_type == b2Body.b2_staticBody || this.m_type == b2Body.b2_kinematicBody) {
      return;
    }
    const center = b2Vec2.Make(0, 0);
    for (let f = this.m_fixtureList; f; f = f.m_next) {
      if (f.m_density == 0.0) {
        continue;
      }
      const massData = f.GetMassData();
      this.m_mass += massData.mass;
      center.x += massData.center.x * massData.mass;
      center.y += massData.center.y * massData.mass;
      this.m_I += massData.I;
    }
    if (this.m_mass > 0.0) {
      this.m_invMass = 1.0 / this.m_mass;
      center.x *= this.m_invMass;
      center.y *= this.m_invMass;
    } else {
      this.m_mass = 1.0;
      this.m_invMass = 1.0;
    }
    if (this.m_I > 0.0 && (this.m_flags & b2Body.e_fixedRotationFlag) == 0) {
      this.m_I -= this.m_mass * (center.x * center.x + center.y * center.y);
      this.m_I *= this.m_inertiaScale;
      b2Settings.b2Assert(this.m_I > 0);
      this.m_invI = 1.0 / this.m_I;
    } else {
      this.m_I = 0.0;
      this.m_invI = 0.0;
    }
    const oldCenter = this.m_sweep.c.Copy();
    this.m_sweep.localCenter.SetV(center);
    this.m_sweep.c0.SetV(b2Math.MulX(this.m_xf, this.m_sweep.localCenter));
    this.m_sweep.c.SetV(this.m_sweep.c0);
    this.m_linearVelocity.x += this.m_angularVelocity * -(this.m_sweep.c.y - oldCenter.y);
    this.m_linearVelocity.y += this.m_angularVelocity * (+(this.m_sweep.c.x - oldCenter.x));
  }
  GetWorldPoint(localPoint) {
    const A = this.m_xf.R;
    const u = new b2Vec2(A.col1.x * localPoint.x + A.col2.x * localPoint.y, A.col1.y * localPoint.x + A.col2.y * localPoint.y);
    u.x += this.m_xf.position.x;
    u.y += this.m_xf.position.y;
    return u;
  }
  GetWorldVector(localVector) {
    return b2Math.MulMV(this.m_xf.R, localVector);
  }
  GetLocalPoint(worldPoint) {
    return b2Math.MulXT(this.m_xf, worldPoint);
  }
  GetLocalVector(worldVector) {
    return b2Math.MulTMV(this.m_xf.R, worldVector);
  }
  GetLinearVelocityFromWorldPoint(worldPoint) {
    return new b2Vec2(this.m_linearVelocity.x - this.m_angularVelocity * (worldPoint.y - this.m_sweep.c.y), this.m_linearVelocity.y + this.m_angularVelocity * (worldPoint.x - this.m_sweep.c.x));
  }
  GetLinearVelocityFromLocalPoint(localPoint) {
    const A = this.m_xf.R;
    const worldPoint = new b2Vec2(A.col1.x * localPoint.x + A.col2.x * localPoint.y, A.col1.y * localPoint.x + A.col2.y * localPoint.y);
    worldPoint.x += this.m_xf.position.x;
    worldPoint.y += this.m_xf.position.y;
    return new b2Vec2(this.m_linearVelocity.x - this.m_angularVelocity * (worldPoint.y - this.m_sweep.c.y), this.m_linearVelocity.y + this.m_angularVelocity * (worldPoint.x - this.m_sweep.c.x));
  }
  GetLinearDamping() {
    return this.m_linearDamping;
  }
  SetLinearDamping(linearDamping = 0) {
    this.m_linearDamping = linearDamping;
  }
  GetAngularDamping() {
    return this.m_angularDamping;
  }
  SetAngularDamping(angularDamping = 0) {
    this.m_angularDamping = angularDamping;
  }
  SetType(type = 0) {
    if (this.m_type == type) {
      return;
    }
    this.m_type = type;
    this.ResetMassData();
    if (this.m_type == b2Body.b2_staticBody) {
      this.m_linearVelocity.SetZero();
      this.m_angularVelocity = 0.0;
    }
    this.SetAwake(true);
    this.m_force.SetZero();
    this.m_torque = 0.0;
    for (let ce = this.m_contactList; ce; ce = ce.next) {
      ce.contact.FlagForFiltering();
    }
  }
  GetType() {
    return this.m_type;
  }
  SetBullet(flag) {
    if (flag) {
      this.m_flags |= b2Body.e_bulletFlag;
    } else {
      this.m_flags &= ~b2Body.e_bulletFlag;
    }
  }
  IsBullet() {
    return (this.m_flags & b2Body.e_bulletFlag) == b2Body.e_bulletFlag;
  }
  SetSleepingAllowed(flag) {
    if (flag) {
      this.m_flags |= b2Body.e_allowSleepFlag;
    } else {
      this.m_flags &= ~b2Body.e_allowSleepFlag;
      this.SetAwake(true);
    }
  }
  SetAwake(flag) {
    if (flag) {
      this.m_flags |= b2Body.e_awakeFlag;
      this.m_sleepTime = 0.0;
    } else {
      this.m_flags &= ~b2Body.e_awakeFlag;
      this.m_sleepTime = 0.0;
      this.m_linearVelocity.SetZero();
      this.m_angularVelocity = 0.0;
      this.m_force.SetZero();
      this.m_torque = 0.0;
    }
  }
  IsAwake() {
    return (this.m_flags & b2Body.e_awakeFlag) == b2Body.e_awakeFlag;
  }
  SetFixedRotation(fixed) {
    if (fixed) {
      this.m_flags |= b2Body.e_fixedRotationFlag;
    } else {
      this.m_flags &= ~b2Body.e_fixedRotationFlag;
    }
    this.ResetMassData();
  }
  IsFixedRotation() {
    return (this.m_flags & b2Body.e_fixedRotationFlag) == b2Body.e_fixedRotationFlag;
  }
  SetActive(flag) {
    if (flag == this.IsActive()) {
      return;
    }
    let broadPhase;
    if (flag) {
      this.m_flags |= b2Body.e_activeFlag;
      broadPhase = this.m_world.m_contactManager.m_broadPhase;
      for (let f = this.m_fixtureList; f; f = f.m_next) {
        f.CreateProxy(broadPhase, this.m_xf);
      }
    } else {
      this.m_flags &= ~b2Body.e_activeFlag;
      broadPhase = this.m_world.m_contactManager.m_broadPhase;
      for (let f = this.m_fixtureList; f; f = f.m_next) {
        f.DestroyProxy(broadPhase);
      }
      let ce = this.m_contactList;
      while (ce) {
        const ce0 = ce;
        ce = ce.next;
        this.m_world.m_contactManager.Destroy(ce0.contact);
      }
      this.m_contactList = null;
    }
  }
  IsActive() {
    return (this.m_flags & b2Body.e_activeFlag) == b2Body.e_activeFlag;
  }
  IsSleepingAllowed() {
    return (this.m_flags & b2Body.e_allowSleepFlag) == b2Body.e_allowSleepFlag;
  }
  GetFixtureList() {
    return this.m_fixtureList;
  }
  GetJointList() {
    return this.m_jointList;
  }
  GetControllerList() {
    return this.m_controllerList;
  }
  GetContactList() {
    return this.m_contactList;
  }
  GetNext() {
    return this.m_next;
  }
  GetUserData() {
    return this.m_userData;
  }
  SetUserData(data) {
    this.m_userData = data;
  }
  GetWorld() {
    return this.m_world;
  }
  SynchronizeFixtures() {
    const xf1 = b2Body.s_xf1;
    const sweep = this.m_sweep;
    const tMat = xf1.R;
    tMat.Set(sweep.a0);
    const tVec = sweep.localCenter;
    xf1.position.x = sweep.c0.x - (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    xf1.position.y = sweep.c0.y - (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
    const broadPhase = this.m_world.m_contactManager.m_broadPhase;
    for (let f = this.m_fixtureList; f; f = f.m_next) {
      f.Synchronize(broadPhase, xf1, this.m_xf);
    }
  }
  SynchronizeTransform() {
    const xf = this.m_xf;
    const sweep = this.m_sweep;
    const tMat = xf.R;
    tMat.Set(sweep.a);
    const tVec = sweep.localCenter;
    xf.position.x = sweep.c.x - (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
    xf.position.y = sweep.c.y - (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
  }
  ShouldCollide(other) {
    if (this.m_type != b2Body.b2_dynamicBody && other.m_type != b2Body.b2_dynamicBody) {
      return false;
    }
    for (let jn = this.m_jointList; jn; jn = jn.next) {
      if (jn.other == other) if (jn.joint.m_collideConnected == false) {
        return false;
      }
    }
    return true;
  }
  Advance(t = 0) {
    const sweep = this.m_sweep;
    sweep.Advance(t);
    sweep.c.SetV(sweep.c0);
    sweep.a = sweep.a0;
    this.SynchronizeTransform();
  }
};

const b2BodyDef = class {
  constructor() {
    this.position = new b2Vec2();
    this.linearVelocity = new b2Vec2();
    this.userData = null;
    this.position.Set(0.0, 0.0);
    this.angle = 0.0;
    this.linearVelocity.Set(0, 0);
    this.angularVelocity = 0.0;
    this.linearDamping = 0.0;
    this.angularDamping = 0.0;
    this.allowSleep = true;
    this.awake = true;
    this.fixedRotation = false;
    this.bullet = false;
    this.type = b2Body.b2_staticBody;
    this.active = true;
    this.inertiaScale = 1.0;
  }
};

const b2ContactFilter = class {
  static b2_defaultFilter = new this();
  ShouldCollide(fixtureA, fixtureB) {
    const filter1 = fixtureA.GetFilterData();
    const filter2 = fixtureB.GetFilterData();
    if (filter1.groupIndex == filter2.groupIndex && filter1.groupIndex != 0) {
      return filter1.groupIndex > 0;
    }
    const collide = (filter1.maskBits & filter2.categoryBits) != 0 && (filter1.categoryBits & filter2.maskBits) != 0;
    return collide;
  }
  RayCollide(userData, fixture) {
    if (!userData) return true;
    return this.ShouldCollide(userData instanceof b2Fixture ? userData : null, fixture);
  }
};

const b2ContactImpulse = class {
  constructor() {
    this.normalImpulses = new Array(b2Settings.b2_maxManifoldPoints).fill(0);
    this.tangentImpulses = new Array(b2Settings.b2_maxManifoldPoints).fill(0);
  }
};

const b2ContactListener = class {
  static b2_defaultListener = new this();
  BeginContact(contact) {}
  EndContact(contact) {}
  PreSolve(contact, oldManifold) {}
  PostSolve(contact, impulse) {}
};

const b2ContactManager = class {
  static s_evalCP = new b2ContactPoint();
  constructor() {
    this.m_world = null;
    this.m_contactCount = 0;
    this.m_contactFilter = b2ContactFilter.b2_defaultFilter;
    this.m_contactListener = b2ContactListener.b2_defaultListener;
    this.m_contactFactory = new b2ContactFactory(this.m_allocator);
    this.m_broadPhase = new b2DynamicTreeBroadPhase();
    this.m_addPairCallback = (proxyUserDataA, proxyUserDataB) => {
      this.AddPair(proxyUserDataA, proxyUserDataB);
    };
  }
  AddPair(proxyUserDataA, proxyUserDataB) {
    let fixtureA = proxyUserDataA;
    let fixtureB = proxyUserDataB;
    let bodyA = fixtureA.m_body;
    let bodyB = fixtureB.m_body;
    if (bodyA == bodyB) return;
    let edge = bodyB.m_contactList;
    while (edge) {
      if (edge.other == bodyA) {
        const fA = edge.contact.m_fixtureA;
        const fB = edge.contact.m_fixtureB;
        if (fA == fixtureA && fB == fixtureB) return;
        if (fA == fixtureB && fB == fixtureA) return;
      }
      edge = edge.next;
    }
    if (bodyB.ShouldCollide(bodyA) == false) {
      return;
    }
    if (this.m_contactFilter.ShouldCollide(fixtureA, fixtureB) == false) {
      return;
    }
    const c = this.m_contactFactory.Create(fixtureA, fixtureB);
    fixtureA = c.m_fixtureA;
    fixtureB = c.m_fixtureB;
    bodyA = fixtureA.m_body;
    bodyB = fixtureB.m_body;
    c.m_prev = null;
    c.m_next = this.m_world.m_contactList;
    if (this.m_world.m_contactList != null) {
      this.m_world.m_contactList.m_prev = c;
    }
    this.m_world.m_contactList = c;
    c.m_nodeA.contact = c;
    c.m_nodeA.other = bodyB;
    c.m_nodeA.prev = null;
    c.m_nodeA.next = bodyA.m_contactList;
    if (bodyA.m_contactList != null) {
      bodyA.m_contactList.prev = c.m_nodeA;
    }
    bodyA.m_contactList = c.m_nodeA;
    c.m_nodeB.contact = c;
    c.m_nodeB.other = bodyA;
    c.m_nodeB.prev = null;
    c.m_nodeB.next = bodyB.m_contactList;
    if (bodyB.m_contactList != null) {
      bodyB.m_contactList.prev = c.m_nodeB;
    }
    bodyB.m_contactList = c.m_nodeB;
    ++this.m_world.m_contactCount;
    return;
  }
  FindNewContacts() {
    this.m_broadPhase.UpdatePairs(this.m_addPairCallback);
  }
  Destroy(c) {
    const bodyA = c.m_fixtureA.m_body;
    const bodyB = c.m_fixtureB.m_body;
    if ((c.m_flags & b2Contact.e_touchingFlag) == b2Contact.e_touchingFlag) {
      this.m_contactListener.EndContact(c);
    }
    if (c.m_prev) {
      c.m_prev.m_next = c.m_next;
    }
    if (c.m_next) {
      c.m_next.m_prev = c.m_prev;
    }
    if (c == this.m_world.m_contactList) {
      this.m_world.m_contactList = c.m_next;
    }
    if (c.m_nodeA.prev) {
      c.m_nodeA.prev.next = c.m_nodeA.next;
    }
    if (c.m_nodeA.next) {
      c.m_nodeA.next.prev = c.m_nodeA.prev;
    }
    if (c.m_nodeA == bodyA.m_contactList) {
      bodyA.m_contactList = c.m_nodeA.next;
    }
    if (c.m_nodeB.prev) {
      c.m_nodeB.prev.next = c.m_nodeB.next;
    }
    if (c.m_nodeB.next) {
      c.m_nodeB.next.prev = c.m_nodeB.prev;
    }
    if (c.m_nodeB == bodyB.m_contactList) {
      bodyB.m_contactList = c.m_nodeB.next;
    }
    this.m_contactFactory.Destroy(c);
    --this.m_world.m_contactCount;
  }
  Collide() {
    let c = this.m_world.m_contactList;
    while (c) {
      const fixtureA = c.m_fixtureA;
      const fixtureB = c.m_fixtureB;
      const bodyA = fixtureA.m_body;
      const bodyB = fixtureB.m_body;
      if ((bodyA.m_flags & b2Body.e_awakeFlag) == 0 && (bodyB.m_flags & b2Body.e_awakeFlag) == 0) {
        c = c.m_next;
        continue;
      }
      if (c.m_flags & b2Contact.e_filterFlag) {
        if (bodyB.ShouldCollide(bodyA) == false) {
          const cNuke = c;
          c = cNuke.m_next;
          this.Destroy(cNuke);
          continue;
        }
        if (this.m_contactFilter.ShouldCollide(fixtureA, fixtureB) == false) {
          const cNuke = c;
          c = cNuke.m_next;
          this.Destroy(cNuke);
          continue;
        }
        c.m_flags &= ~b2Contact.e_filterFlag;
      }
      const proxyA = fixtureA.m_proxy;
      const proxyB = fixtureB.m_proxy;
      const overlap = this.m_broadPhase.TestOverlap(proxyA, proxyB);
      if (overlap == false) {
        const cNuke = c;
        c = cNuke.m_next;
        this.Destroy(cNuke);
        continue;
      }
      c.Update(this.m_contactListener);
      c = c.m_next;
    }
  }
};

const b2DebugDraw = class {
  static e_shapeBit = 0x0001;
  static e_jointBit = 0x0002;
  static e_aabbBit = 0x0004;
  static e_pairBit = 0x0008;
  static e_centerOfMassBit = 0x0010;
  static e_controllerBit = 0x0020;
  constructor() {
    this.m_drawScale = 1.0;
    this.m_lineThickness = 1.0;
    this.m_alpha = 1.0;
    this.m_fillAlpha = 1.0;
    this.m_xformScale = 1.0;
    this.m_sprite = {
      graphics: {
        clear: () => {
          this.m_ctx.clearRect(0, 0, this.m_ctx.canvas.width, this.m_ctx.canvas.height);
        }
      }
    };
    this.m_drawFlags = 0;
  }
  _color(color, alpha) {
    return "rgba(" + ((color & 0xFF0000) >> 16) + "," + ((color & 0xFF00) >> 8) + "," + (color & 0xFF) + "," + alpha + ")";
  }
  SetFlags(flags = 0) {
    this.m_drawFlags = flags;
  }
  GetFlags() {
    return this.m_drawFlags;
  }
  AppendFlags(flags = 0) {
    this.m_drawFlags |= flags;
  }
  ClearFlags(flags = 0) {
    this.m_drawFlags &= ~flags;
  }
  SetSprite(sprite) {
    this.m_ctx = sprite;
  }
  GetSprite() {
    return this.m_ctx;
  }
  SetDrawScale(drawScale = 0) {
    this.m_drawScale = drawScale;
  }
  GetDrawScale() {
    return this.m_drawScale;
  }
  SetLineThickness(lineThickness = 0) {
    this.m_lineThickness = lineThickness;
    this.m_ctx.strokeWidth = lineThickness;
  }
  GetLineThickness() {
    return this.m_lineThickness;
  }
  SetAlpha(alpha = 0) {
    this.m_alpha = alpha;
  }
  GetAlpha() {
    return this.m_alpha;
  }
  SetFillAlpha(alpha = 0) {
    this.m_fillAlpha = alpha;
  }
  GetFillAlpha() {
    return this.m_fillAlpha;
  }
  SetXFormScale(xformScale = 0) {
    this.m_xformScale = xformScale;
  }
  GetXFormScale() {
    return this.m_xformScale;
  }
  DrawPolygon(vertices, vertexCount, color) {
    if (!vertexCount) return;
    const s = this.m_ctx;
    const drawScale = this.m_drawScale;
    s.beginPath();
    s.strokeStyle = this._color(color.color, this.m_alpha);
    s.moveTo(vertices[0].x * drawScale, vertices[0].y * drawScale);
    for (let i = 1; i < vertexCount; i++) {
      s.lineTo(vertices[i].x * drawScale, vertices[i].y * drawScale);
    }
    s.lineTo(vertices[0].x * drawScale, vertices[0].y * drawScale);
    s.closePath();
    s.stroke();
  }
  DrawSolidPolygon(vertices, vertexCount, color) {
    if (!vertexCount) return;
    const s = this.m_ctx;
    const drawScale = this.m_drawScale;
    s.beginPath();
    s.strokeStyle = this._color(color.color, this.m_alpha);
    s.fillStyle = this._color(color.color, this.m_fillAlpha);
    s.moveTo(vertices[0].x * drawScale, vertices[0].y * drawScale);
    for (let i = 1; i < vertexCount; i++) {
      s.lineTo(vertices[i].x * drawScale, vertices[i].y * drawScale);
    }
    s.lineTo(vertices[0].x * drawScale, vertices[0].y * drawScale);
    s.closePath();
    s.fill();
    s.stroke();
  }
  DrawCircle(center, radius, color) {
    if (!radius) return;
    const s = this.m_ctx;
    const drawScale = this.m_drawScale;
    s.beginPath();
    s.strokeStyle = this._color(color.color, this.m_alpha);
    s.arc(center.x * drawScale, center.y * drawScale, radius * drawScale, 0, Math.PI * 2, true);
    s.closePath();
    s.stroke();
  }
  DrawSolidCircle(center, radius, axis, color) {
    if (!radius) return;
    const s = this.m_ctx,
      drawScale = this.m_drawScale,
      cx = center.x * drawScale,
      cy = center.y * drawScale;
    s.moveTo(0, 0);
    s.beginPath();
    s.strokeStyle = this._color(color.color, this.m_alpha);
    s.fillStyle = this._color(color.color, this.m_fillAlpha);
    s.arc(cx, cy, radius * drawScale, 0, Math.PI * 2, true);
    s.moveTo(cx, cy);
    s.lineTo((center.x + axis.x * radius) * drawScale, (center.y + axis.y * radius) * drawScale);
    s.closePath();
    s.fill();
    s.stroke();
  }
  DrawSegment(p1, p2, color) {
    const s = this.m_ctx,
      drawScale = this.m_drawScale;
    s.strokeStyle = this._color(color.color, this.m_alpha);
    s.beginPath();
    s.moveTo(p1.x * drawScale, p1.y * drawScale);
    s.lineTo(p2.x * drawScale, p2.y * drawScale);
    s.closePath();
    s.stroke();
  }
  DrawTransform(xf) {
    const s = this.m_ctx,
      drawScale = this.m_drawScale;
    s.beginPath();
    s.strokeStyle = this._color(0xff0000, this.m_alpha);
    s.moveTo(xf.position.x * drawScale, xf.position.y * drawScale);
    s.lineTo((xf.position.x + this.m_xformScale * xf.R.col1.x) * drawScale, (xf.position.y + this.m_xformScale * xf.R.col1.y) * drawScale);
    s.strokeStyle = this._color(0xff00, this.m_alpha);
    s.moveTo(xf.position.x * drawScale, xf.position.y * drawScale);
    s.lineTo((xf.position.x + this.m_xformScale * xf.R.col2.x) * drawScale, (xf.position.y + this.m_xformScale * xf.R.col2.y) * drawScale);
    s.closePath();
    s.stroke();
  }
};

const b2DestructionListener = class {
  SayGoodbyeJoint(joint) {}
  SayGoodbyeFixture(fixture) {}
};

const b2FilterData = class {
  constructor() {
    this.categoryBits = 0x0001;
    this.maskBits = 0xFFFF;
    this.groupIndex = 0;
  }
  Copy() {
    const copy = new b2FilterData();
    copy.categoryBits = this.categoryBits;
    copy.maskBits = this.maskBits;
    copy.groupIndex = this.groupIndex;
    return copy;
  }
};

const b2Fixture = class {
  static s_aabb1 = new b2AABB();
  static s_aabb2 = new b2AABB();
  static s_displacement = new b2Vec2();
  constructor() {
    this.m_filter = new b2FilterData();
    this.m_aabb = new b2AABB();
    this.m_userData = null;
    this.m_body = null;
    this.m_next = null;
    this.m_shape = null;
    this.m_density = 0.0;
    this.m_friction = 0.0;
    this.m_restitution = 0.0;
  }
  GetType() {
    return this.m_shape.GetType();
  }
  GetShape() {
    return this.m_shape;
  }
  SetSensor(sensor) {
    if (this.m_isSensor == sensor) return;
    this.m_isSensor = sensor;
    if (this.m_body == null) return;
    let edge = this.m_body.GetContactList();
    while (edge) {
      const contact = edge.contact;
      const fixtureA = contact.GetFixtureA();
      const fixtureB = contact.GetFixtureB();
      if (fixtureA == this || fixtureB == this) contact.SetSensor(fixtureA.IsSensor() || fixtureB.IsSensor());
      edge = edge.next;
    }
  }
  IsSensor() {
    return this.m_isSensor;
  }
  SetFilterData(filter) {
    this.m_filter = filter.Copy();
    if (!this.m_body) return;
    let edge = this.m_body.GetContactList();
    while (edge) {
      const contact = edge.contact;
      const fixtureA = contact.GetFixtureA();
      const fixtureB = contact.GetFixtureB();
      if (fixtureA == this || fixtureB == this) contact.FlagForFiltering();
      edge = edge.next;
    }
  }
  GetFilterData() {
    return this.m_filter.Copy();
  }
  GetBody() {
    return this.m_body;
  }
  GetNext() {
    return this.m_next;
  }
  GetUserData() {
    return this.m_userData;
  }
  SetUserData(data) {
    this.m_userData = data;
  }
  TestPoint(p) {
    return this.m_shape.TestPoint(this.m_body.GetTransform(), p);
  }
  RayCast(output, input) {
    return this.m_shape.RayCast(output, input, this.m_body.GetTransform());
  }
  GetMassData(massData = null) {
    if (massData == null) {
      massData = new b2MassData();
    }
    this.m_shape.ComputeMass(massData, this.m_density);
    return massData;
  }
  SetDensity(density = 0) {
    this.m_density = density;
  }
  GetDensity() {
    return this.m_density;
  }
  GetFriction() {
    return this.m_friction;
  }
  SetFriction(friction = 0) {
    this.m_friction = friction;
  }
  GetRestitution() {
    return this.m_restitution;
  }
  SetRestitution(restitution = 0) {
    this.m_restitution = restitution;
  }
  GetAABB() {
    return this.m_aabb;
  }
  Create(body, xf, def) {
    this.m_userData = def.userData;
    this.m_friction = def.friction;
    this.m_restitution = def.restitution;
    this.m_body = body;
    this.m_next = null;
    this.m_filter = def.filter.Copy();
    this.m_isSensor = def.isSensor;
    this.m_shape = def.shape.Copy();
    this.m_density = def.density;
  }
  Destroy() {
    this.m_shape = null;
  }
  CreateProxy(broadPhase, xf) {
    this.m_shape.ComputeAABB(this.m_aabb, xf);
    this.m_proxy = broadPhase.CreateProxy(this.m_aabb, this);
  }
  DestroyProxy(broadPhase) {
    if (this.m_proxy == null) {
      return;
    }
    broadPhase.DestroyProxy(this.m_proxy);
    this.m_proxy = null;
  }
  Synchronize(broadPhase, transform1, transform2) {
    if (!this.m_proxy) return;
    const aabb1 = b2Fixture.s_aabb1;
    const aabb2 = b2Fixture.s_aabb2;
    this.m_shape.ComputeAABB(aabb1, transform1);
    this.m_shape.ComputeAABB(aabb2, transform2);
    this.m_aabb.Combine(aabb1, aabb2);
    const displacement = b2Fixture.s_displacement;
    displacement.Set(transform2.position.x - transform1.position.x, transform2.position.y - transform1.position.y);
    broadPhase.MoveProxy(this.m_proxy, this.m_aabb, displacement);
  }
};

const b2FixtureDef = class {
  constructor() {
    this.filter = new b2FilterData();
    this.shape = null;
    this.userData = null;
    this.friction = 0.2;
    this.restitution = 0.0;
    this.density = 0.0;
    this.filter.categoryBits = 0x0001;
    this.filter.maskBits = 0xFFFF;
    this.filter.groupIndex = 0;
    this.isSensor = false;
  }
};

const b2Island = class {
  static s_impulse = new b2ContactImpulse();
  constructor() {
    this.m_bodies = [];
    this.m_contacts = [];
    this.m_joints = [];
  }
  Initialize(bodyCapacity = 0, contactCapacity = 0, jointCapacity = 0, allocator, listener, contactSolver) {
    this.m_bodyCapacity = bodyCapacity;
    this.m_contactCapacity = contactCapacity;
    this.m_jointCapacity = jointCapacity;
    this.m_bodyCount = 0;
    this.m_contactCount = 0;
    this.m_jointCount = 0;
    this.m_allocator = allocator;
    this.m_listener = listener;
    this.m_contactSolver = contactSolver;
    for (let i = this.m_bodies.length; i < bodyCapacity; i++) {
      this.m_bodies[i] = null;
    }
    for (let i = this.m_contacts.length; i < contactCapacity; i++) {
      this.m_contacts[i] = null;
    }
    for (let i = this.m_joints.length; i < jointCapacity; i++) {
      this.m_joints[i] = null;
    }
  }
  Clear() {
    this.m_bodyCount = 0;
    this.m_contactCount = 0;
    this.m_jointCount = 0;
  }
  Solve(step, gravity, allowSleep) {
    let b;
    let joint;
    const dt = step.dt;
    const bodies = this.m_bodies;
    const bodyCount = this.m_bodyCount;
    let linearVelocity;
    let sweep;
    for (let i = 0; i < bodyCount; ++i) {
      b = bodies[i];
      if (b.m_type != b2Body.b2_dynamicBody) continue;
      linearVelocity = b.m_linearVelocity;
      linearVelocity.x += dt * (gravity.x + b.m_invMass * b.m_force.x);
      linearVelocity.y += dt * (gravity.y + b.m_invMass * b.m_force.y);
      b.m_angularVelocity += dt * b.m_invI * b.m_torque;
      linearVelocity.Multiply(b2Math.Clamp(1.0 - dt * b.m_linearDamping, 0.0, 1.0));
      b.m_angularVelocity *= b2Math.Clamp(1.0 - dt * b.m_angularDamping, 0.0, 1.0);
    }
    this.m_contactSolver.Initialize(step, this.m_contacts, this.m_contactCount, this.m_allocator);
    const contactSolver = this.m_contactSolver;
    contactSolver.InitVelocityConstraints(step);
    for (let i = 0; i < this.m_jointCount; ++i) {
      joint = this.m_joints[i];
      joint.InitVelocityConstraints(step);
    }
    for (let i = 0; i < step.velocityIterations; ++i) {
      for (let j = 0; j < this.m_jointCount; ++j) {
        joint = this.m_joints[j];
        joint.SolveVelocityConstraints(step);
      }
      contactSolver.SolveVelocityConstraints();
    }
    for (let i = 0; i < this.m_jointCount; ++i) {
      joint = this.m_joints[i];
      joint.FinalizeVelocityConstraints();
    }
    contactSolver.FinalizeVelocityConstraints();
    for (let i = 0; i < bodyCount; ++i) {
      b = bodies[i];
      if (b.m_type == b2Body.b2_staticBody) continue;
      linearVelocity = b.m_linearVelocity;
      const translationX = dt * linearVelocity.x;
      const translationY = dt * linearVelocity.y;
      if ((translationX * translationX + translationY * translationY) > b2Settings.b2_maxTranslationSquared) {
        linearVelocity.Normalize();
        linearVelocity.x *= b2Settings.b2_maxTranslation * step.inv_dt;
        linearVelocity.y *= b2Settings.b2_maxTranslation * step.inv_dt;
      }
      const rotation = dt * b.m_angularVelocity;
      if (rotation * rotation > b2Settings.b2_maxRotationSquared) {
        if (b.m_angularVelocity < 0.0) {
          b.m_angularVelocity = -b2Settings.b2_maxRotation * step.inv_dt;
        } else {
          b.m_angularVelocity = b2Settings.b2_maxRotation * step.inv_dt;
        }
      }
      sweep = b.m_sweep;
      sweep.c0.SetV(sweep.c);
      sweep.a0 = sweep.a;
      sweep.c.x += dt * linearVelocity.x;
      sweep.c.y += dt * linearVelocity.y;
      sweep.a += dt * b.m_angularVelocity;
      b.SynchronizeTransform();
    }
    for (let i = 0; i < step.positionIterations; ++i) {
      const contactsOkay = contactSolver.SolvePositionConstraints(b2Settings.b2_contactBaumgarte);
      let jointsOkay = true;
      for (let j = 0; j < this.m_jointCount; ++j) {
        joint = this.m_joints[j];
        const jointOkay = joint.SolvePositionConstraints(b2Settings.b2_contactBaumgarte);
        jointsOkay = jointsOkay && jointOkay;
      }
      if (contactsOkay && jointsOkay) {
        break;
      }
    }
    this.Report(contactSolver.m_constraints);
    if (allowSleep) {
      let minSleepTime = Number.MAX_VALUE;
      const linTolSqr = b2Settings.b2_linearSleepTolerance * b2Settings.b2_linearSleepTolerance;
      const angTolSqr = b2Settings.b2_angularSleepTolerance * b2Settings.b2_angularSleepTolerance;
      for (let i = 0; i < bodyCount; ++i) {
        b = bodies[i];
        if (b.m_type == b2Body.b2_staticBody) {
          continue;
        }
        if ((b.m_flags & b2Body.e_allowSleepFlag) == 0) {
          b.m_sleepTime = 0.0;
          minSleepTime = 0.0;
        }
        if ((b.m_flags & b2Body.e_allowSleepFlag) == 0 || b.m_angularVelocity * b.m_angularVelocity > angTolSqr || b2Math.Dot(b.m_linearVelocity, b.m_linearVelocity) > linTolSqr) {
          b.m_sleepTime = 0.0;
          minSleepTime = 0.0;
        } else {
          b.m_sleepTime += dt;
          minSleepTime = b2Math.Min(minSleepTime, b.m_sleepTime);
        }
      }
      if (minSleepTime >= b2Settings.b2_timeToSleep) {
        for (let i = 0; i < bodyCount; ++i) {
          b = bodies[i];
          b.SetAwake(false);
        }
      }
    }
  }
  SolveTOI(subStep) {
    let b;
    const dt = subStep.dt;
    const bodies = this.m_bodies;
    const bodyCount = this.m_bodyCount;
    let linearVelocity;
    let sweep;
    this.m_contactSolver.Initialize(subStep, this.m_contacts, this.m_contactCount, this.m_allocator);
    const contactSolver = this.m_contactSolver;
    contactSolver.InitVelocityConstraints(subStep);
    for (let i = 0; i < this.m_jointCount; ++i) {
      this.m_joints[i].InitVelocityConstraints(subStep);
    }
    for (let i = 0; i < subStep.velocityIterations; ++i) {
      contactSolver.SolveVelocityConstraints();
      for (let j = 0; j < this.m_jointCount; ++j) {
        this.m_joints[j].SolveVelocityConstraints(subStep);
      }
    }
    for (let i = 0; i < bodyCount; ++i) {
      b = bodies[i];
      if (b.m_type == b2Body.b2_staticBody) continue;
      linearVelocity = b.m_linearVelocity;
      const translationX = dt * linearVelocity.x;
      const translationY = dt * linearVelocity.y;
      if ((translationX * translationX + translationY * translationY) > b2Settings.b2_maxTranslationSquared) {
        linearVelocity.Normalize();
        linearVelocity.x *= b2Settings.b2_maxTranslation * subStep.inv_dt;
        linearVelocity.y *= b2Settings.b2_maxTranslation * subStep.inv_dt;
      }
      const rotation = dt * b.m_angularVelocity;
      if (rotation * rotation > b2Settings.b2_maxRotationSquared) {
        if (b.m_angularVelocity < 0.0) {
          b.m_angularVelocity = -b2Settings.b2_maxRotation * subStep.inv_dt;
        } else {
          b.m_angularVelocity = b2Settings.b2_maxRotation * subStep.inv_dt;
        }
      }
      sweep = b.m_sweep;
      sweep.c0.SetV(sweep.c);
      sweep.a0 = sweep.a;
      sweep.c.x += dt * linearVelocity.x;
      sweep.c.y += dt * linearVelocity.y;
      sweep.a += dt * b.m_angularVelocity;
      b.SynchronizeTransform();
    }
    const k_toiBaumgarte = 0.75;
    for (let i = 0; i < subStep.positionIterations; ++i) {
      const contactsOkay = contactSolver.SolvePositionConstraints(k_toiBaumgarte);
      let jointsOkay = true;
      for (let j = 0; j < this.m_jointCount; ++j) {
        const jointOkay = this.m_joints[j].SolvePositionConstraints(b2Settings.b2_contactBaumgarte);
        jointsOkay = jointsOkay && jointOkay;
      }
      if (contactsOkay && jointsOkay) {
        break;
      }
    }
    this.Report(contactSolver.m_constraints);
  }
  Report(constraints) {
    if (this.m_listener == null) {
      return;
    }
    for (let i = 0; i < this.m_contactCount; ++i) {
      const c = this.m_contacts[i];
      const cc = constraints[i];
      for (let j = 0; j < cc.pointCount; ++j) {
        b2Island.s_impulse.normalImpulses[j] = cc.points[j].normalImpulse;
        b2Island.s_impulse.tangentImpulses[j] = cc.points[j].tangentImpulse;
      }
      this.m_listener.PostSolve(c, b2Island.s_impulse);
    }
  }
  AddBody(body) {
    body.m_islandIndex = this.m_bodyCount;
    this.m_bodies[this.m_bodyCount++] = body;
  }
  AddContact(contact) {
    this.m_contacts[this.m_contactCount++] = contact;
  }
  AddJoint(joint) {
    this.m_joints[this.m_jointCount++] = joint;
  }
};

const b2TimeStep = class {
  constructor() {
    this.dt = 0.0;
    this.inv_dt = 0.0;
    this.dtRatio = 0.0;
    this.velocityIterations = 0;
    this.positionIterations = 0;
    this.warmStarting = false;
  }
  Set(step) {
    this.dt = step.dt;
    this.inv_dt = step.inv_dt;
    this.positionIterations = step.positionIterations;
    this.velocityIterations = step.velocityIterations;
    this.warmStarting = step.warmStarting;
  }
};

const b2World = class {
  static s_queue = [];
  static e_newFixture = 0x0001;
  static e_locked = 0x0002;
  static s_timestep2 = new b2TimeStep();
  static s_xf = new b2Transform();
  static s_backupA = new b2Sweep();
  static s_backupB = new b2Sweep();
  static s_timestep = new b2TimeStep();
  static s_jointColor = new b2Color(0.5, 0.8, 0.8);
  constructor(gravity, doSleep) {
    this.s_stack = [];
    this.m_contactManager = new b2ContactManager();
    this.m_contactSolver = new b2ContactSolver();
    this.m_island = new b2Island();
    this.m_destructionListener = null;
    this.m_debugDraw = null;
    this.m_bodyList = null;
    this.m_contactList = null;
    this.m_jointList = null;
    this.m_controllerList = null;
    this.m_bodyCount = 0;
    this.m_contactCount = 0;
    this.m_jointCount = 0;
    this.m_controllerCount = 0;
    this.m_warmStarting = true;
    this.m_continuousPhysics = true;
    this.m_allowSleep = doSleep;
    this.m_gravity = gravity;
    this.m_inv_dt0 = 0.0;
    this.m_stepRemainder = 0.0;
    this.m_contactManager.m_world = this;
    const bd = new b2BodyDef();
    this.m_groundBody = this.CreateBody(bd);
  }
  SetDestructionListener(listener) {
    this.m_destructionListener = listener;
  }
  SetContactFilter(filter) {
    this.m_contactManager.m_contactFilter = filter;
  }
  SetContactListener(listener) {
    this.m_contactManager.m_contactListener = listener;
  }
  SetDebugDraw(debugDraw) {
    this.m_debugDraw = debugDraw;
  }
  SetBroadPhase(broadPhase) {
    const oldBroadPhase = this.m_contactManager.m_broadPhase;
    this.m_contactManager.m_broadPhase = broadPhase;
    for (let b = this.m_bodyList; b; b = b.m_next) {
      for (let f = b.m_fixtureList; f; f = f.m_next) {
        f.m_proxy = broadPhase.CreateProxy(oldBroadPhase.GetFatAABB(f.m_proxy), f);
      }
    }
  }
  Validate() {
    this.m_contactManager.m_broadPhase.Validate();
  }
  GetProxyCount() {
    return this.m_contactManager.m_broadPhase.GetProxyCount();
  }
  CreateBody(def) {
    if (this.IsLocked() == true) {
      return null;
    }
    const b = new b2Body(def, this);
    b.m_prev = null;
    b.m_next = this.m_bodyList;
    if (this.m_bodyList) {
      this.m_bodyList.m_prev = b;
    }
    this.m_bodyList = b;
    ++this.m_bodyCount;
    return b;
  }
  DestroyBody(b) {
    if (this.IsLocked() == true) {
      return;
    }
    let jn = b.m_jointList;
    while (jn) {
      const jn0 = jn;
      jn = jn.next;
      if (this.m_destructionListener) {
        this.m_destructionListener.SayGoodbyeJoint(jn0.joint);
      }
      this.DestroyJoint(jn0.joint);
    }
    let coe = b.m_controllerList;
    while (coe) {
      const coe0 = coe;
      coe = coe.nextController;
      coe0.controller.RemoveBody(b);
    }
    let ce = b.m_contactList;
    while (ce) {
      const ce0 = ce;
      ce = ce.next;
      this.m_contactManager.Destroy(ce0.contact);
    }
    b.m_contactList = null;
    let f = b.m_fixtureList;
    while (f) {
      const f0 = f;
      f = f.m_next;
      if (this.m_destructionListener) {
        this.m_destructionListener.SayGoodbyeFixture(f0);
      }
      f0.DestroyProxy(this.m_contactManager.m_broadPhase);
      f0.Destroy();
    }
    b.m_fixtureList = null;
    b.m_fixtureCount = 0;
    if (b.m_prev) {
      b.m_prev.m_next = b.m_next;
    }
    if (b.m_next) {
      b.m_next.m_prev = b.m_prev;
    }
    if (b == this.m_bodyList) {
      this.m_bodyList = b.m_next;
    }
    --this.m_bodyCount;
  }
  CreateJoint(def) {
    const j = b2Joint.Create(def, null);
    j.m_prev = null;
    j.m_next = this.m_jointList;
    if (this.m_jointList) {
      this.m_jointList.m_prev = j;
    }
    this.m_jointList = j;
    ++this.m_jointCount;
    j.m_edgeA.joint = j;
    j.m_edgeA.other = j.m_bodyB;
    j.m_edgeA.prev = null;
    j.m_edgeA.next = j.m_bodyA.m_jointList;
    if (j.m_bodyA.m_jointList) j.m_bodyA.m_jointList.prev = j.m_edgeA;
    j.m_bodyA.m_jointList = j.m_edgeA;
    j.m_edgeB.joint = j;
    j.m_edgeB.other = j.m_bodyA;
    j.m_edgeB.prev = null;
    j.m_edgeB.next = j.m_bodyB.m_jointList;
    if (j.m_bodyB.m_jointList) j.m_bodyB.m_jointList.prev = j.m_edgeB;
    j.m_bodyB.m_jointList = j.m_edgeB;
    const bodyA = def.bodyA;
    const bodyB = def.bodyB;
    if (def.collideConnected == false) {
      let edge = bodyB.GetContactList();
      while (edge) {
        if (edge.other == bodyA) {
          edge.contact.FlagForFiltering();
        }
        edge = edge.next;
      }
    }
    return j;
  }
  DestroyJoint(j) {
    const collideConnected = j.m_collideConnected;
    if (j.m_prev) {
      j.m_prev.m_next = j.m_next;
    }
    if (j.m_next) {
      j.m_next.m_prev = j.m_prev;
    }
    if (j == this.m_jointList) {
      this.m_jointList = j.m_next;
    }
    const bodyA = j.m_bodyA;
    const bodyB = j.m_bodyB;
    bodyA.SetAwake(true);
    bodyB.SetAwake(true);
    if (j.m_edgeA.prev) {
      j.m_edgeA.prev.next = j.m_edgeA.next;
    }
    if (j.m_edgeA.next) {
      j.m_edgeA.next.prev = j.m_edgeA.prev;
    }
    if (j.m_edgeA == bodyA.m_jointList) {
      bodyA.m_jointList = j.m_edgeA.next;
    }
    j.m_edgeA.prev = null;
    j.m_edgeA.next = null;
    if (j.m_edgeB.prev) {
      j.m_edgeB.prev.next = j.m_edgeB.next;
    }
    if (j.m_edgeB.next) {
      j.m_edgeB.next.prev = j.m_edgeB.prev;
    }
    if (j.m_edgeB == bodyB.m_jointList) {
      bodyB.m_jointList = j.m_edgeB.next;
    }
    j.m_edgeB.prev = null;
    j.m_edgeB.next = null;
    b2Joint.Destroy(j, null);
    --this.m_jointCount;
    if (collideConnected == false) {
      let edge = bodyB.GetContactList();
      while (edge) {
        if (edge.other == bodyA) {
          edge.contact.FlagForFiltering();
        }
        edge = edge.next;
      }
    }
  }
  AddController(c) {
    c.m_next = this.m_controllerList;
    c.m_prev = null;
    this.m_controllerList = c;
    c.m_world = this;
    this.m_controllerCount++;
    return c;
  }
  RemoveController(c) {
    if (c.m_prev) c.m_prev.m_next = c.m_next;
    if (c.m_next) c.m_next.m_prev = c.m_prev;
    if (this.m_controllerList == c) this.m_controllerList = c.m_next;
    this.m_controllerCount--;
  }
  CreateController(controller) {
    if (controller.m_world != this) throw new Error("Controller can only be a member of one world");
    controller.m_next = this.m_controllerList;
    controller.m_prev = null;
    if (this.m_controllerList) this.m_controllerList.m_prev = controller;
    this.m_controllerList = controller;
    ++this.m_controllerCount;
    controller.m_world = this;
    return controller;
  }
  DestroyController(controller) {
    controller.Clear();
    if (controller.m_next) controller.m_next.m_prev = controller.m_prev;
    if (controller.m_prev) controller.m_prev.m_next = controller.m_next;
    if (controller == this.m_controllerList) this.m_controllerList = controller.m_next;
    --this.m_controllerCount;
  }
  SetWarmStarting(flag) {
    this.m_warmStarting = flag;
  }
  SetContinuousPhysics(flag) {
    this.m_continuousPhysics = flag;
  }
  GetBodyCount() {
    return this.m_bodyCount;
  }
  GetJointCount() {
    return this.m_jointCount;
  }
  GetContactCount() {
    return this.m_contactCount;
  }
  SetGravity(gravity) {
    this.m_gravity = gravity;
  }
  GetGravity() {
    return this.m_gravity;
  }
  GetGroundBody() {
    return this.m_groundBody;
  }
  Step(dt = 0, velocityIterations = 0, positionIterations = 0) {
    if (this.m_flags & b2World.e_newFixture) {
      this.m_contactManager.FindNewContacts();
      this.m_flags &= ~b2World.e_newFixture;
    }
    this.m_flags |= b2World.e_locked;
    const step = b2World.s_timestep2;
    step.dt = dt;
    step.velocityIterations = velocityIterations;
    step.positionIterations = positionIterations;
    if (dt > 0.0) {
      step.inv_dt = 1.0 / dt;
    } else {
      step.inv_dt = 0.0;
    }
    step.dtRatio = this.m_inv_dt0 * dt;
    step.warmStarting = this.m_warmStarting;
    this.m_contactManager.Collide();
    if (step.dt > 0.0) {
      this.Solve(step);
    }
    if (this.m_continuousPhysics && step.dt > 0.0) {
      this.SolveTOI(step);
    }
    if (step.dt > 0.0) {
      this.m_inv_dt0 = step.inv_dt;
    }
    this.m_flags &= ~b2World.e_locked;
  }
  StepFixed(dt = 0, timeStep = 0, velocityIterations = 0, positionIterations = 0) {
    b2Settings.b2Assert(timeStep > 0.0);
    let stepCount = 0;
    if (dt > 0.0 && b2Math.IsValid(dt)) {
      this.m_stepRemainder += dt;
    }
    while (this.m_stepRemainder >= timeStep) {
      this.Step(timeStep, velocityIterations, positionIterations);
      this.m_stepRemainder -= timeStep;
      ++stepCount;
    }
    return stepCount;
  }
  GetStepRemainder() {
    return this.m_stepRemainder;
  }
  ClearForces() {
    for (let body = this.m_bodyList; body; body = body.m_next) {
      body.m_force.SetZero();
      body.m_torque = 0.0;
    }
  }
  DrawDebugData() {
    if (this.m_debugDraw == null) {
      return;
    }
    this.m_debugDraw.m_sprite.graphics.clear();
    const flags = this.m_debugDraw.GetFlags();
    let s;
    let bp;
    let xf;
    let vs;
    const color = new b2Color(0, 0, 0);
    if (flags & b2DebugDraw.e_shapeBit) {
      for (let b = this.m_bodyList; b; b = b.m_next) {
        xf = b.m_xf;
        for (let f = b.GetFixtureList(); f; f = f.m_next) {
          s = f.GetShape();
          if (b.IsActive() == false) {
            color.Set(0.5, 0.5, 0.3);
            this.DrawShape(s, xf, color);
          } else if (b.GetType() == b2Body.b2_staticBody) {
            color.Set(0.5, 0.9, 0.5);
            this.DrawShape(s, xf, color);
          } else if (b.GetType() == b2Body.b2_kinematicBody) {
            color.Set(0.5, 0.5, 0.9);
            this.DrawShape(s, xf, color);
          } else if (b.IsAwake() == false) {
            color.Set(0.6, 0.6, 0.6);
            this.DrawShape(s, xf, color);
          } else {
            color.Set(0.9, 0.7, 0.7);
            this.DrawShape(s, xf, color);
          }
        }
      }
    }
    if (flags & b2DebugDraw.e_jointBit) {
      for (let j = this.m_jointList; j; j = j.m_next) {
        this.DrawJoint(j);
      }
    }
    if (flags & b2DebugDraw.e_controllerBit) {
      for (let c = this.m_controllerList; c; c = c.m_next) {
        c.Draw(this.m_debugDraw);
      }
    }
    if (flags & b2DebugDraw.e_pairBit) {
      color.Set(0.3, 0.9, 0.9);
      for (let contact = this.m_contactManager.m_contactList; contact; contact = contact.GetNext()) {
        const fixtureA = contact.GetFixtureA();
        const fixtureB = contact.GetFixtureB();
        const cA = fixtureA.GetAABB().GetCenter();
        const cB = fixtureB.GetAABB().GetCenter();
        this.m_debugDraw.DrawSegment(cA, cB, color);
      }
    }
    if (flags & b2DebugDraw.e_aabbBit) {
      bp = this.m_contactManager.m_broadPhase;
      vs = [new b2Vec2(), new b2Vec2(), new b2Vec2(), new b2Vec2()];
      for (let b = this.m_bodyList; b; b = b.GetNext()) {
        if (b.IsActive() == false) {
          continue;
        }
        for (let f = b.GetFixtureList(); f; f = f.GetNext()) {
          const aabb = bp.GetFatAABB(f.m_proxy);
          vs[0].Set(aabb.lowerBound.x, aabb.lowerBound.y);
          vs[1].Set(aabb.upperBound.x, aabb.lowerBound.y);
          vs[2].Set(aabb.upperBound.x, aabb.upperBound.y);
          vs[3].Set(aabb.lowerBound.x, aabb.upperBound.y);
          this.m_debugDraw.DrawPolygon(vs, 4, color);
        }
      }
    }
    if (flags & b2DebugDraw.e_centerOfMassBit) {
      for (let b = this.m_bodyList; b; b = b.m_next) {
        xf = b2World.s_xf;
        xf.R = b.m_xf.R;
        xf.position = b.GetWorldCenter();
        this.m_debugDraw.DrawTransform(xf);
      }
    }
  }
  QueryAABB(callback, aabb) {
    const broadPhase = this.m_contactManager.m_broadPhase;
    const WorldQueryWrapper = (proxy) => {
      return callback(broadPhase.GetUserData(proxy));
    };
    broadPhase.Query(WorldQueryWrapper, aabb);
  }
  QueryShape(callback, shape, transform = null) {
    if (transform == null) {
      transform = new b2Transform();
      transform.SetIdentity();
    }
    const broadPhase = this.m_contactManager.m_broadPhase;
    const WorldQueryWrapper = (proxy) => {
      const fixture = broadPhase.GetUserData(proxy);
      if (b2Shape.TestOverlap(shape, transform, fixture.GetShape(), fixture.GetBody().GetTransform())) return callback(fixture);
      return true;
    };
    const aabb = new b2AABB();
    shape.ComputeAABB(aabb, transform);
    broadPhase.Query(WorldQueryWrapper, aabb);
  }
  QueryPoint(callback, p) {
    const broadPhase = this.m_contactManager.m_broadPhase;
    const WorldQueryWrapper = (proxy) => {
      const fixture = broadPhase.GetUserData(proxy);
      if (fixture.TestPoint(p)) return callback(fixture);
      return true;
    };
    const aabb = new b2AABB();
    aabb.lowerBound.Set(p.x - b2Settings.b2_linearSlop, p.y - b2Settings.b2_linearSlop);
    aabb.upperBound.Set(p.x + b2Settings.b2_linearSlop, p.y + b2Settings.b2_linearSlop);
    broadPhase.Query(WorldQueryWrapper, aabb);
  }
  RayCast(callback, point1, point2) {
    const broadPhase = this.m_contactManager.m_broadPhase;
    const output = new b2RayCastOutput();
    const RayCastWrapper = (input, proxy) => {
      const userData = broadPhase.GetUserData(proxy);
      const fixture = userData;
      if (fixture.RayCast(output, input)) {
        const fraction = output.fraction;
        const point = new b2Vec2((1.0 - fraction) * point1.x + fraction * point2.x, (1.0 - fraction) * point1.y + fraction * point2.y);
        return callback(fixture, point, output.normal, fraction);
      }
      return input.maxFraction;
    };
    const input = new b2RayCastInput(point1, point2);
    broadPhase.RayCast(RayCastWrapper, input);
  }
  RayCastOne(point1, point2) {
    const broadPhase = this.m_contactManager.m_broadPhase;
    const output = new b2RayCastOutput();
    let result;
    const RayCastOneWrapper = (input, proxy) => {
      const userData = broadPhase.GetUserData(proxy);
      const fixture = userData;
      if (fixture.RayCast(output, input)) {
        result = fixture;
        return output.fraction;
      }
      return input.maxFraction;
    };
    const input = new b2RayCastInput(point1, point2);
    broadPhase.RayCast(RayCastOneWrapper, input);
    return result;
  }
  RayCastAll(point1, point2) {
    const broadPhase = this.m_contactManager.m_broadPhase;
    const output = new b2RayCastOutput();
    const result = [];
    const RayCastAllWrapper = (input, proxy) => {
      const userData = broadPhase.GetUserData(proxy);
      const fixture = userData;
      if (fixture.RayCast(output, input)) {
        result[result.length] = fixture;
        return 1;
      }
      return input.maxFraction;
    };
    const input = new b2RayCastInput(point1, point2);
    broadPhase.RayCast(RayCastAllWrapper, input);
    return result;
  }
  GetBodyList() {
    return this.m_bodyList;
  }
  GetJointList() {
    return this.m_jointList;
  }
  GetContactList() {
    return this.m_contactList;
  }
  IsLocked() {
    return (this.m_flags & b2World.e_locked) > 0;
  }
  Solve(step) {
    let b;
    for (let controller = this.m_controllerList; controller; controller = controller.m_next) {
      controller.Step(step);
    }
    const island = this.m_island;
    island.Initialize(this.m_bodyCount, this.m_contactCount, this.m_jointCount, null, this.m_contactManager.m_contactListener, this.m_contactSolver);
    for (b = this.m_bodyList; b; b = b.m_next) {
      b.m_flags &= ~b2Body.e_islandFlag;
    }
    for (let c = this.m_contactList; c; c = c.m_next) {
      c.m_flags &= ~b2Contact.e_islandFlag;
    }
    for (let j = this.m_jointList; j; j = j.m_next) {
      j.m_islandFlag = false;
    }
    const stack = this.s_stack;
    for (let seed = this.m_bodyList; seed; seed = seed.m_next) {
      if (seed.m_flags & b2Body.e_islandFlag) {
        continue;
      }
      if ((seed.m_flags & b2Body.e_awakeFlag) == 0 || (seed.m_flags & b2Body.e_activeFlag) == 0) {
        continue;
      }
      if (seed.m_type == b2Body.b2_staticBody) {
        continue;
      }
      island.Clear();
      let stackCount = 0;
      stack[stackCount++] = seed;
      seed.m_flags |= b2Body.e_islandFlag;
      while (stackCount > 0) {
        b = stack[--stackCount];
        island.AddBody(b);
        if ((b.m_flags & b2Body.e_awakeFlag) == 0) {
          b.SetAwake(true);
        }
        if (b.m_type == b2Body.b2_staticBody) {
          continue;
        }
        for (let ce = b.m_contactList; ce; ce = ce.next) {
          const contact = ce.contact;
          const contactFlags = contact.m_flags;
          if (contactFlags & b2Contact.e_islandFlag) {
            continue;
          }
          if ((contactFlags & b2Contact.e_sensorFlag) == b2Contact.e_sensorFlag ||(contactFlags & b2Contact.e_enabledFlag) == 0 || (contactFlags & b2Contact.e_touchingFlag) == 0) {
            continue;
          }
          island.AddContact(contact);
          contact.m_flags |= b2Contact.e_islandFlag;
          const other = ce.other;
          if (other.m_flags & b2Body.e_islandFlag) {
            continue;
          }
          stack[stackCount++] = other;
          other.m_flags |= b2Body.e_islandFlag;
        }
        for (let jn = b.m_jointList; jn; jn = jn.next) {
          if (jn.joint.m_islandFlag == true) {
            continue;
          }
          const other = jn.other;
          if ((other.m_flags & b2Body.e_activeFlag) == 0) {
            continue;
          }
          island.AddJoint(jn.joint);
          jn.joint.m_islandFlag = true;
          if (other.m_flags & b2Body.e_islandFlag) {
            continue;
          }
          stack[stackCount++] = other;
          other.m_flags |= b2Body.e_islandFlag;
        }
      }
      island.Solve(step, this.m_gravity, this.m_allowSleep);
      for (let i = 0; i < island.m_bodyCount; ++i) {
        b = island.m_bodies[i];
        if (b.m_type == b2Body.b2_staticBody) {
          b.m_flags &= ~b2Body.e_islandFlag;
        }
      }
    }
    for (let i = 0; i < stack.length; ++i) {
      if (!stack[i]) break;
      stack[i] = null;
    }
    for (b = this.m_bodyList; b; b = b.m_next) {
      if ((b.m_flags & b2Body.e_awakeFlag) == 0 || (b.m_flags & b2Body.e_activeFlag) == 0) {
        continue;
      }
      if (b.m_type == b2Body.b2_staticBody) {
        continue;
      }
      b.SynchronizeFixtures();
    }
    this.m_contactManager.FindNewContacts();
  }
  SolveTOI(step) {
    let b;
    let fA;
    let fB;
    let bA;
    let bB;
    let j;
    let sweepA;
    let sweepB;
    const island = this.m_island;
    island.Initialize(this.m_bodyCount, b2Settings.b2_maxTOIContactsPerIsland, b2Settings.b2_maxTOIJointsPerIsland, null, this.m_contactManager.m_contactListener, this.m_contactSolver);
    const queue = b2World.s_queue;
    for (b = this.m_bodyList; b; b = b.m_next) {
      b.m_flags &= ~b2Body.e_islandFlag;
      b.m_sweep.t0 = 0.0;
    }
    let c;
    for (c = this.m_contactList; c; c = c.m_next) {
      c.m_flags &= ~ (b2Contact.e_toiFlag | b2Contact.e_islandFlag);
    }
    for (j = this.m_jointList; j; j = j.m_next) {
      j.m_islandFlag = false;
    }
    for (;;) {
      let minContact = null;
      let minTOI = 1.0;
      for (c = this.m_contactList; c; c = c.m_next) {
        if ((c.m_flags & b2Contact.e_sensorFlag) == b2Contact.e_sensorFlag || (c.m_flags & b2Contact.e_enabledFlag) == 0 || (c.m_flags & b2Contact.e_continuousFlag) == 0) {
          continue;
        }
        let toi;
        if (c.m_flags & b2Contact.e_toiFlag) {
          toi = c.m_toi;
        } else {
          bA = c.m_fixtureA.m_body;
          bB = c.m_fixtureB.m_body;
          if ((bA.m_type != b2Body.b2_dynamicBody || (bA.m_flags & b2Body.e_awakeFlag) == 0) && (bB.m_type != b2Body.b2_dynamicBody || (bB.m_flags & b2Body.e_awakeFlag) == 0)) {
            continue;
          }
          sweepA = bA.m_sweep;
          sweepB = bB.m_sweep;
          let t0 = sweepA.t0;
          if (sweepA.t0 < sweepB.t0) {
            t0 = sweepB.t0;
            sweepA.Advance(t0);
          } else if (sweepB.t0 < sweepA.t0) {
            t0 = sweepA.t0;
            sweepB.Advance(t0);
          }
          toi = c.ComputeTOI(sweepA, sweepB);
          b2Settings.b2Assert(0.0 <= toi && toi <= 1.0);
          if (toi > 0.0 && toi < 1.0) {
            toi = (1.0 - toi) * t0 + toi;
            if (toi > 1) toi = 1;
          }
          c.m_toi = toi;
          c.m_flags |= b2Contact.e_toiFlag;
        }
        if (b2Settings.b2_epsilon < toi && toi < minTOI) {
          minContact = c;
          minTOI = toi;
        }
      }
      if (minContact == null || 1.0 - 100.0 * b2Settings.b2_epsilon < minTOI) {
        break;
      }
      fA = minContact.m_fixtureA;
      fB = minContact.m_fixtureB;
      bA = fA.m_body;
      bB = fB.m_body;
      b2World.s_backupA.Set(bA.m_sweep);
      b2World.s_backupB.Set(bB.m_sweep);
      bA.Advance(minTOI);
      bB.Advance(minTOI);
      minContact.Update(this.m_contactManager.m_contactListener);
      minContact.m_flags &= ~b2Contact.e_toiFlag;
      if ((minContact.m_flags & b2Contact.e_sensorFlag) == b2Contact.e_sensorFlag || (minContact.m_flags & b2Contact.e_enabledFlag) == 0) {
        bA.m_sweep.Set(b2World.s_backupA);
        bB.m_sweep.Set(b2World.s_backupB);
        bA.SynchronizeTransform();
        bB.SynchronizeTransform();
        continue;
      }
      if ((minContact.m_flags & b2Contact.e_touchingFlag) == 0) {
        continue;
      }
      let seed = bA;
      if (seed.m_type != b2Body.b2_dynamicBody) {
        seed = bB;
      }
      island.Clear();
      let queueStart = 0;
      let queueSize = 0;
      queue[queueStart + queueSize++] = seed;
      seed.m_flags |= b2Body.e_islandFlag;
      while (queueSize > 0) {
        b = queue[queueStart++];
        --queueSize;
        island.AddBody(b);
        if ((b.m_flags & b2Body.e_awakeFlag) == 0) {
          b.SetAwake(true);
        }
        if (b.m_type != b2Body.b2_dynamicBody) {
          continue;
        }
        for (let cEdge = b.m_contactList; cEdge; cEdge = cEdge.next) {
          if (island.m_contactCount == island.m_contactCapacity) {
            break;
          }
          const contact = cEdge.contact;
          const contactFlags = contact.m_flags;
          if (contactFlags & b2Contact.e_islandFlag) {
            continue;
          }
          if ((contactFlags & b2Contact.e_sensorFlag) == b2Contact.e_sensorFlag || (contactFlags & b2Contact.e_enabledFlag) == 0 || (contactFlags & b2Contact.e_touchingFlag) == 0) {
            continue;
          }
          island.AddContact(contact);
          contact.m_flags |= b2Contact.e_islandFlag;
          const other = cEdge.other;
          if (other.m_flags & b2Body.e_islandFlag) {
            continue;
          }
          if (other.m_type != b2Body.b2_staticBody) {
            other.Advance(minTOI);
            other.SetAwake(true);
          }
          queue[queueStart + queueSize] = other;
          ++queueSize;
          other.m_flags |= b2Body.e_islandFlag;
        }
        for (let jEdge = b.m_jointList; jEdge; jEdge = jEdge.next) {
          if (island.m_jointCount == island.m_jointCapacity) continue;
          if (jEdge.joint.m_islandFlag == true) continue;
          const other = jEdge.other;
          if ((other.m_flags & b2Body.e_activeFlag) == 0) {
            continue;
          }
          island.AddJoint(jEdge.joint);
          jEdge.joint.m_islandFlag = true;
          if (other.m_flags & b2Body.e_islandFlag) continue;
          if (other.m_type != b2Body.b2_staticBody) {
            other.Advance(minTOI);
            other.SetAwake(true);
          }
          queue[queueStart + queueSize] = other;
          ++queueSize;
          other.m_flags |= b2Body.e_islandFlag;
        }
      }
      const subStep = b2World.s_timestep;
      subStep.warmStarting = false;
      subStep.dt = (1.0 - minTOI) * step.dt;
      subStep.inv_dt = 1.0 / subStep.dt;
      subStep.dtRatio = 0.0;
      subStep.velocityIterations = step.velocityIterations;
      subStep.positionIterations = step.positionIterations;
      island.SolveTOI(subStep);
      for (let i = 0; i < island.m_bodyCount; ++i) {
        b = island.m_bodies[i];
        b.m_flags &= ~b2Body.e_islandFlag;
        if ((b.m_flags & b2Body.e_awakeFlag) == 0) {
          continue;
        }
        if (b.m_type != b2Body.b2_dynamicBody) {
          continue;
        }
        b.SynchronizeFixtures();
        for (let cEdge = b.m_contactList; cEdge; cEdge = cEdge.next) {
          cEdge.contact.m_flags &= ~b2Contact.e_toiFlag;
        }
      }
      for (let i = 0; i < island.m_contactCount; ++i) {
        c = island.m_contacts[i];
        c.m_flags &= ~ (b2Contact.e_toiFlag | b2Contact.e_islandFlag);
      }
      for (let i = 0; i < island.m_jointCount; ++i) {
        j = island.m_joints[i];
        j.m_islandFlag = false;
      }
      this.m_contactManager.FindNewContacts();
    }
  }
  DrawJoint(joint) {
    const b1 = joint.GetBodyA();
    const b2 = joint.GetBodyB();
    const xf1 = b1.m_xf;
    const xf2 = b2.m_xf;
    const x1 = xf1.position;
    const x2 = xf2.position;
    const p1 = joint.GetAnchorA();
    const p2 = joint.GetAnchorB();
    const color = b2World.s_jointColor;
    switch (joint.m_type) {
    case b2Joint.e_distanceJoint:
      this.m_debugDraw.DrawSegment(p1, p2, color);
      break;
    case b2Joint.e_pulleyJoint: {
      const pulley = joint;
      const s1 = pulley.GetGroundAnchorA();
      const s2 = pulley.GetGroundAnchorB();
      this.m_debugDraw.DrawSegment(s1, p1, color);
      this.m_debugDraw.DrawSegment(s2, p2, color);
      this.m_debugDraw.DrawSegment(s1, s2, color);
    }
      break;
    case b2Joint.e_mouseJoint:
      this.m_debugDraw.DrawSegment(p1, p2, color);
      break;
    default:
      if (b1 != this.m_groundBody) this.m_debugDraw.DrawSegment(x1, p1, color);
      this.m_debugDraw.DrawSegment(p1, p2, color);
      if (b2 != this.m_groundBody) this.m_debugDraw.DrawSegment(x2, p2, color);
    }
  }
  DrawShape(shape, xf, color) {
    switch (shape.m_type) {
    case b2Shape.e_circleShape: {
      const circle = shape;
      const center = b2Math.MulX(xf, circle.m_p);
      const radius = circle.m_radius;
      const axis = xf.R.col1;
      this.m_debugDraw.DrawSolidCircle(center, radius, axis, color);
    }
      break;
    case b2Shape.e_polygonShape: {
      const poly = shape;
      const vertexCount = poly.GetVertexCount();
      const localVertices = poly.GetVertices();
      const vertices = new Array(vertexCount);
      for (let i = 0; i < vertexCount; ++i) {
        vertices[i] = b2Math.MulX(xf, localVertices[i]);
      }
      this.m_debugDraw.DrawSolidPolygon(vertices, vertexCount, color);
    }
      break;
    case b2Shape.e_edgeShape: {
      const edge = shape;
      this.m_debugDraw.DrawSegment(b2Math.MulX(xf, edge.GetVertex1()), b2Math.MulX(xf, edge.GetVertex2()), color);
    }
      break;
    }
  }
};

const b2Contact = class {
  static e_sensorFlag = 0x0001;
  static e_continuousFlag = 0x0002;
  static e_islandFlag = 0x0004;
  static e_toiFlag = 0x0008;
  static e_touchingFlag = 0x0010;
  static e_enabledFlag = 0x0020;
  static e_filterFlag = 0x0040;
  static s_input = new b2TOIInput();
  constructor() {
    this.m_nodeA = new b2ContactEdge();
    this.m_nodeB = new b2ContactEdge();
    this.m_manifold = new b2Manifold();
    this.m_oldManifold = new b2Manifold();
    this.m_flags = 0;
    this.m_fixtureA = null;
    this.m_fixtureB = null;
    this.m_prev = null;
    this.m_next = null;
    this.m_toi = 0.0;
  }
  GetManifold() {
    return this.m_manifold;
  }
  GetWorldManifold(worldManifold) {
    const bodyA = this.m_fixtureA.GetBody();
    const bodyB = this.m_fixtureB.GetBody();
    const shapeA = this.m_fixtureA.GetShape();
    const shapeB = this.m_fixtureB.GetShape();
    worldManifold.Initialize(this.m_manifold, bodyA.GetTransform(), shapeA.m_radius, bodyB.GetTransform(), shapeB.m_radius);
  }
  IsTouching() {
    return (this.m_flags & b2Contact.e_touchingFlag) == b2Contact.e_touchingFlag;
  }
  IsContinuous() {
    return (this.m_flags & b2Contact.e_continuousFlag) == b2Contact.e_continuousFlag;
  }
  SetSensor(sensor) {
    if (sensor) {
      this.m_flags |= b2Contact.e_sensorFlag;
    } else {
      this.m_flags &= ~b2Contact.e_sensorFlag;
    }
  }
  IsSensor() {
    return (this.m_flags & b2Contact.e_sensorFlag) == b2Contact.e_sensorFlag;
  }
  SetEnabled(flag) {
    if (flag) {
      this.m_flags |= b2Contact.e_enabledFlag;
    } else {
      this.m_flags &= ~b2Contact.e_enabledFlag;
    }
  }
  IsEnabled() {
    return (this.m_flags & b2Contact.e_enabledFlag) == b2Contact.e_enabledFlag;
  }
  GetNext() {
    return this.m_next;
  }
  GetFixtureA() {
    return this.m_fixtureA;
  }
  GetFixtureB() {
    return this.m_fixtureB;
  }
  FlagForFiltering() {
    this.m_flags |= b2Contact.e_filterFlag;
  }
  Reset(fixtureA = null, fixtureB = null) {
    this.m_flags = b2Contact.e_enabledFlag;
    if (!fixtureA || !fixtureB) {
      this.m_fixtureA = null;
      this.m_fixtureB = null;
      return;
    }
    if (fixtureA.IsSensor() || fixtureB.IsSensor()) {
      this.m_flags |= b2Contact.e_sensorFlag;
    }
    const bodyA = fixtureA.GetBody();
    const bodyB = fixtureB.GetBody();
    if (bodyA.GetType() != b2Body.b2_dynamicBody || bodyA.IsBullet() || bodyB.GetType() != b2Body.b2_dynamicBody || bodyB.IsBullet()) {
      this.m_flags |= b2Contact.e_continuousFlag;
    }
    this.m_fixtureA = fixtureA;
    this.m_fixtureB = fixtureB;
    this.m_manifold.m_pointCount = 0;
    this.m_prev = null;
    this.m_next = null;
    this.m_nodeA.contact = null;
    this.m_nodeA.prev = null;
    this.m_nodeA.next = null;
    this.m_nodeA.other = null;
    this.m_nodeB.contact = null;
    this.m_nodeB.prev = null;
    this.m_nodeB.next = null;
    this.m_nodeB.other = null;
  }
  Update(listener) {
    const tManifold = this.m_oldManifold;
    this.m_oldManifold = this.m_manifold;
    this.m_manifold = tManifold;
    this.m_flags |= b2Contact.e_enabledFlag;
    let touching = false;
    const wasTouching = (this.m_flags & b2Contact.e_touchingFlag) == b2Contact.e_touchingFlag;
    const fixtureA = this.m_fixtureA;
    const fixtureB = this.m_fixtureB;
    const bodyA = fixtureA.m_body;
    const bodyB = fixtureB.m_body;
    const aabbOverlap = fixtureA.m_aabb.TestOverlap(fixtureB.m_aabb);
    if (this.m_flags & b2Contact.e_sensorFlag) {
      if (aabbOverlap) {
        touching = b2Shape.TestOverlap(fixtureA.m_shape, bodyA.m_xf, fixtureB.m_shape, bodyB.m_xf);
      }
      this.m_manifold.m_pointCount = 0;
    } else {
      if (bodyA.m_type != b2Body.b2_dynamicBody || (bodyA.m_flags & b2Body.e_bulletFlag) == b2Body.e_bulletFlag || bodyB.m_type != b2Body.b2_dynamicBody || (bodyB.m_flags & b2Body.e_bulletFlag) == b2Body.e_bulletFlag) {
        this.m_flags |= b2Contact.e_continuousFlag;
      } else {
        this.m_flags &= ~b2Contact.e_continuousFlag;
      }
      if (aabbOverlap) {
        this.Evaluate();
        const manifold = this.m_manifold;
        const pointCount = manifold.m_pointCount;
        const oldPoints = this.m_oldManifold.m_points;
        const oldPointCount = this.m_oldManifold.m_pointCount;
        touching = pointCount > 0;
        for (let i = 0; i < pointCount; ++i) {
          const mp2 = manifold.m_points[i];
          mp2.m_normalImpulse = 0.0;
          mp2.m_tangentImpulse = 0.0;
          const key2 = mp2.m_id._key;
          for (let j = 0; j < oldPointCount; ++j) {
            const mp1 = oldPoints[j];
            if (mp1.m_id._key == key2) {
              mp2.m_normalImpulse = mp1.m_normalImpulse;
              mp2.m_tangentImpulse = mp1.m_tangentImpulse;
              break;
            }
          }
        }
      } else {
        this.m_manifold.m_pointCount = 0;
      }
      if (touching != wasTouching) {
        bodyA.SetAwake(true);
        bodyB.SetAwake(true);
      }
    }
    if (touching) {
      this.m_flags |= b2Contact.e_touchingFlag;
    } else {
      this.m_flags &= ~b2Contact.e_touchingFlag;
    }
    if (wasTouching == false && touching == true) {
      listener.BeginContact(this);
    }
    if (wasTouching == true && touching == false) {
      listener.EndContact(this);
    }
    if ((this.m_flags & b2Contact.e_sensorFlag) == 0) {
      listener.PreSolve(this, this.m_oldManifold);
    }
  }
  Evaluate() {}
  ComputeTOI(sweepA, sweepB) {
    const input = b2Contact.s_input;
    input.proxyA.Set(this.m_fixtureA.m_shape);
    input.proxyB.Set(this.m_fixtureB.m_shape);
    input.sweepA = sweepA;
    input.sweepB = sweepB;
    input.tolerance = b2Settings.b2_linearSlop;
    return b2TimeOfImpact.TimeOfImpact(input);
  }
};

const b2CircleContact = class extends b2Contact {
  static Create(allocator) {
    return new b2CircleContact();
  }
  static Destroy(contact, allocator) {}
  Reset(fixtureA, fixtureB) {
    super.Reset(fixtureA, fixtureB);
  }
  Evaluate() {
    const fixtureA = this.m_fixtureA;
    const fixtureB = this.m_fixtureB;
    b2Collision.CollideCircles(this.m_manifold, fixtureA.m_shape, fixtureA.m_body.m_xf, fixtureB.m_shape, fixtureB.m_body.m_xf);
  }
};

const b2ContactConstraint = class {
  constructor() {
    this.localPlaneNormal = new b2Vec2();
    this.localPoint = new b2Vec2();
    this.normal = new b2Vec2();
    this.normalMass = new b2Mat22();
    this.K = new b2Mat22();
    this.bodyA = null;
    this.bodyB = null;
    this.manifold = null;
    this.pointCount = 0;
    this.friction = 0.0;
    this.restitution = 0.0;
    this.radius = 0.0;
    this.type = 0;
    this.points = new Array(b2Settings.b2_maxManifoldPoints);
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.points[i] = new b2ContactConstraintPoint();
    }
  }
};

const b2ContactConstraintPoint = class {
  constructor() {
    this.localPoint = new b2Vec2();
    this.rA = new b2Vec2();
    this.rB = new b2Vec2();
    this.normalImpulse = 0.0;
    this.tangentImpulse = 0.0;
    this.normalMass = 0.0;
    this.equalizedMass = 0.0;
    this.tangentMass = 0.0;
    this.velocityBias = 0.0;
  }
};

const b2ContactEdge = class {
};

const b2ContactFactory = class {
  constructor(allocator) {
    this.m_allocator = allocator;
    this.InitializeRegisters();
  }
  AddType(createFcn, destroyFcn, type1 = 0, type2 = 0) {
    this.m_registers[type1][type2].createFcn = createFcn;
    this.m_registers[type1][type2].destroyFcn = destroyFcn;
    this.m_registers[type1][type2].primary = true;
    if (type1 != type2) {
      this.m_registers[type2][type1].createFcn = createFcn;
      this.m_registers[type2][type1].destroyFcn = destroyFcn;
      this.m_registers[type2][type1].primary = false;
    }
  }
  InitializeRegisters() {
    this.m_registers = new Array(b2Shape.e_shapeTypeCount);
    for (let i = 0; i < b2Shape.e_shapeTypeCount; i++) {
      this.m_registers[i] = new Array(b2Shape.e_shapeTypeCount);
      for (let j = 0; j < b2Shape.e_shapeTypeCount; j++) {
        this.m_registers[i][j] = new b2ContactRegister();
      }
    }
    this.AddType(b2CircleContact.Create, b2CircleContact.Destroy, b2Shape.e_circleShape, b2Shape.e_circleShape);
    this.AddType(b2PolyAndCircleContact.Create, b2PolyAndCircleContact.Destroy, b2Shape.e_polygonShape, b2Shape.e_circleShape);
    this.AddType(b2PolygonContact.Create, b2PolygonContact.Destroy, b2Shape.e_polygonShape, b2Shape.e_polygonShape);
    this.AddType(b2EdgeAndCircleContact.Create, b2EdgeAndCircleContact.Destroy, b2Shape.e_edgeShape, b2Shape.e_circleShape);
    this.AddType(b2PolyAndEdgeContact.Create, b2PolyAndEdgeContact.Destroy, b2Shape.e_polygonShape, b2Shape.e_edgeShape);
  }
  Create(fixtureA, fixtureB) {
    const type1 = fixtureA.m_shape.m_type;
    const type2 = fixtureB.m_shape.m_type;
    const reg = this.m_registers[type1][type2];
    let c;
    if (reg.pool) {
      c = reg.pool;
      reg.pool = c.m_next;
      reg.poolCount--;
      c.Reset(fixtureA, fixtureB);
      return c;
    }
    const createFcn = reg.createFcn;
    if (createFcn != null) {
      if (reg.primary) {
        c = createFcn(this.m_allocator);
        c.Reset(fixtureA, fixtureB);
        return c;
      } else {
        c = createFcn(this.m_allocator);
        c.Reset(fixtureB, fixtureA);
        return c;
      }
    } else {
      return null;
    }
  }
  Destroy(contact) {
    if (contact.m_manifold.m_pointCount > 0) {
      contact.m_fixtureA.m_body.SetAwake(true);
      contact.m_fixtureB.m_body.SetAwake(true);
    }
    const type1 = contact.m_fixtureA.m_shape.m_type;
    const type2 = contact.m_fixtureB.m_shape.m_type;
    const reg = this.m_registers[type1][type2];
    reg.poolCount++;
    contact.m_next = reg.pool;
    reg.pool = contact;
    const destroyFcn = reg.destroyFcn;
    destroyFcn(contact, this.m_allocator);
  }
};

const b2ContactRegister = class {
};

const b2ContactResult = class {
  constructor() {
    this.position = new b2Vec2();
    this.normal = new b2Vec2();
    this.id = new b2ContactID();
  }
};

const b2PositionSolverManifold = class {
  static circlePointA = new b2Vec2();
  static circlePointB = new b2Vec2();
  constructor() {
    this.m_normal = new b2Vec2();
    this.m_separations = new Array(b2Settings.b2_maxManifoldPoints).fill(0);
    this.m_points = new Array(b2Settings.b2_maxManifoldPoints);
    for (let i = 0; i < b2Settings.b2_maxManifoldPoints; i++) {
      this.m_points[i] = new b2Vec2();
    }
  }
  Initialize(cc) {
    b2Settings.b2Assert(cc.pointCount > 0);
    let clipPointX;
    let clipPointY;
    let tMat;
    let tVec;
    let planePointX;
    let planePointY;
    let normalX;
    let normalY;
    const xfA = cc.bodyA.m_xf;
    const xfB = cc.bodyB.m_xf;
    const normal = this.m_normal;
    const manifoldPoints = this.m_points;
    const separations = this.m_separations;
    const pointCount = cc.pointCount;
    const radius = cc.radius;
    const constraintPoints = cc.points;
    let point;
    switch (cc.type) {
    case b2Manifold.e_circles: {
      tMat = xfA.R;
      tVec = cc.localPoint;
      const pointAX = xfA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      const pointAY = xfA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tMat = xfB.R;
      tVec = constraintPoints[0].localPoint;
      const pointBX = xfB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      const pointBY = xfB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      const dX = pointBX - pointAX;
      const dY = pointBY - pointAY;
      const d2 = dX * dX + dY * dY;
      if (d2 > b2Settings.b2_epsilon * b2Settings.b2_epsilon) {
        const d = Math.sqrt(d2);
        normalX = dX / d;
        normalY = dY / d;
      } else {
        normalX = 1.0;
        normalY = 0.0;
      }
      point = manifoldPoints[0];
      point.x = 0.5 * (pointAX + pointBX);
      point.y = 0.5 * (pointAY + pointBY);
      separations[0] = dX * normalX + dY * normalY - radius;
      normal.x = normalX;
      normal.y = normalY;
    }
      break;
    case b2Manifold.e_faceA:
      tMat = xfA.R;
      tVec = cc.localPlaneNormal;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = cc.localPoint;
      planePointX = xfA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      planePointY = xfA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tMat = xfB.R;
      for (let i = 0; i < pointCount; ++i) {
        tVec = constraintPoints[i].localPoint;
        clipPointX = xfB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
        clipPointY = xfB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
        separations[i] = (clipPointX - planePointX) * normalX + (clipPointY - planePointY) * normalY - radius;
        point = manifoldPoints[i];
        point.x = clipPointX;
        point.y = clipPointY;
      }
      normal.x = normalX;
      normal.y = normalY;
      break;
    case b2Manifold.e_faceB:
      tMat = xfB.R;
      tVec = cc.localPlaneNormal;
      normalX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      normalY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tVec = cc.localPoint;
      planePointX = xfB.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
      planePointY = xfB.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
      tMat = xfA.R;
      for (let i = 0; i < pointCount; ++i) {
        tVec = constraintPoints[i].localPoint;
        clipPointX = xfA.position.x + (tMat.col1.x * tVec.x + tMat.col2.x * tVec.y);
        clipPointY = xfA.position.y + (tMat.col1.y * tVec.x + tMat.col2.y * tVec.y);
        separations[i] = (clipPointX - planePointX) * normalX + (clipPointY - planePointY) * normalY - radius;
        point = manifoldPoints[i];
        point.x = clipPointX;
        point.y = clipPointY;
      }
      normal.x = normalX * -1;
      normal.y = normalY * -1;
      break;
    }
  }
};

const b2ContactSolver = class {
  static s_worldManifold = new b2WorldManifold();
  static s_psm = new b2PositionSolverManifold();
  constructor() {
    this.m_step = new b2TimeStep();
    this.m_constraints = [];
    this.m_allocator = null;
    this.m_constraintCount = 0;
  }
  Initialize(step, contacts, contactCount = 0, allocator) {
    let contact;
    this.m_step.Set(step);
    this.m_allocator = allocator;
    this.m_constraintCount = contactCount;
    while (this.m_constraints.length < this.m_constraintCount) {
      this.m_constraints[this.m_constraints.length] = new b2ContactConstraint();
    }
    for (let i = 0; i < contactCount; ++i) {
      contact = contacts[i];
      const fixtureA = contact.m_fixtureA;
      const fixtureB = contact.m_fixtureB;
      const shapeA = fixtureA.m_shape;
      const shapeB = fixtureB.m_shape;
      const radiusA = shapeA.m_radius;
      const radiusB = shapeB.m_radius;
      const bodyA = fixtureA.m_body;
      const bodyB = fixtureB.m_body;
      const manifold = contact.GetManifold();
      const friction = b2Settings.b2MixFriction(fixtureA.GetFriction(), fixtureB.GetFriction());
      const restitution = b2Settings.b2MixRestitution(fixtureA.GetRestitution(), fixtureB.GetRestitution());
      const vAX = bodyA.m_linearVelocity.x;
      const vAY = bodyA.m_linearVelocity.y;
      const vBX = bodyB.m_linearVelocity.x;
      const vBY = bodyB.m_linearVelocity.y;
      const wA = bodyA.m_angularVelocity;
      const wB = bodyB.m_angularVelocity;
      b2Settings.b2Assert(manifold.m_pointCount > 0);
      b2ContactSolver.s_worldManifold.Initialize(manifold, bodyA.m_xf, radiusA, bodyB.m_xf, radiusB);
      const normalX = b2ContactSolver.s_worldManifold.m_normal.x;
      const normalY = b2ContactSolver.s_worldManifold.m_normal.y;
      const cc = this.m_constraints[i];
      cc.bodyA = bodyA;
      cc.bodyB = bodyB;
      cc.manifold = manifold;
      cc.normal.x = normalX;
      cc.normal.y = normalY;
      cc.pointCount = manifold.m_pointCount;
      cc.friction = friction;
      cc.restitution = restitution;
      cc.localPlaneNormal.x = manifold.m_localPlaneNormal.x;
      cc.localPlaneNormal.y = manifold.m_localPlaneNormal.y;
      cc.localPoint.x = manifold.m_localPoint.x;
      cc.localPoint.y = manifold.m_localPoint.y;
      cc.radius = radiusA + radiusB;
      cc.type = manifold.m_type;
      for (let k = 0; k < cc.pointCount; ++k) {
        const cp = manifold.m_points[k];
        const ccp = cc.points[k];
        ccp.normalImpulse = cp.m_normalImpulse;
        ccp.tangentImpulse = cp.m_tangentImpulse;
        ccp.localPoint.SetV(cp.m_localPoint);
        const rAX = ccp.rA.x = b2ContactSolver.s_worldManifold.m_points[k].x - bodyA.m_sweep.c.x;
        const rAY = ccp.rA.y = b2ContactSolver.s_worldManifold.m_points[k].y - bodyA.m_sweep.c.y;
        const rBX = ccp.rB.x = b2ContactSolver.s_worldManifold.m_points[k].x - bodyB.m_sweep.c.x;
        const rBY = ccp.rB.y = b2ContactSolver.s_worldManifold.m_points[k].y - bodyB.m_sweep.c.y;
        let rnA = rAX * normalY - rAY * normalX;
        let rnB = rBX * normalY - rBY * normalX;
        rnA *= rnA;
        rnB *= rnB;
        const kNormal = bodyA.m_invMass + bodyB.m_invMass + bodyA.m_invI * rnA + bodyB.m_invI * rnB;
        ccp.normalMass = 1.0 / kNormal;
        let kEqualized = bodyA.m_mass * bodyA.m_invMass + bodyB.m_mass * bodyB.m_invMass;
        kEqualized += bodyA.m_mass * bodyA.m_invI * rnA + bodyB.m_mass * bodyB.m_invI * rnB;
        ccp.equalizedMass = 1.0 / kEqualized;
        const tangentX = normalY;
        const tangentY = -normalX;
        let rtA = rAX * tangentY - rAY * tangentX;
        let rtB = rBX * tangentY - rBY * tangentX;
        rtA *= rtA;
        rtB *= rtB;
        const kTangent = bodyA.m_invMass + bodyB.m_invMass + bodyA.m_invI * rtA + bodyB.m_invI * rtB;
        ccp.tangentMass = 1.0 / kTangent;
        ccp.velocityBias = 0.0;
        const tX = vBX + (-wB * rBY) - vAX - (-wA * rAY);
        const tY = vBY + (wB * rBX) - vAY - (wA * rAX);
        const vRel = cc.normal.x * tX + cc.normal.y * tY;
        if (vRel < -b2Settings.b2_velocityThreshold) {
          ccp.velocityBias += (-cc.restitution * vRel);
        }
      }
      if (cc.pointCount == 2) {
        const ccp1 = cc.points[0];
        const ccp2 = cc.points[1];
        const invMassA = bodyA.m_invMass;
        const invIA = bodyA.m_invI;
        const invMassB = bodyB.m_invMass;
        const invIB = bodyB.m_invI;
        const rn1A = ccp1.rA.x * normalY - ccp1.rA.y * normalX;
        const rn1B = ccp1.rB.x * normalY - ccp1.rB.y * normalX;
        const rn2A = ccp2.rA.x * normalY - ccp2.rA.y * normalX;
        const rn2B = ccp2.rB.x * normalY - ccp2.rB.y * normalX;
        const k11 = invMassA + invMassB + invIA * rn1A * rn1A + invIB * rn1B * rn1B;
        const k22 = invMassA + invMassB + invIA * rn2A * rn2A + invIB * rn2B * rn2B;
        const k12 = invMassA + invMassB + invIA * rn1A * rn2A + invIB * rn1B * rn2B;
        const k_maxConditionNumber = 100.0;
        if (k11 * k11 < k_maxConditionNumber * (k11 * k22 - k12 * k12)) {
          cc.K.col1.Set(k11, k12);
          cc.K.col2.Set(k12, k22);
          cc.K.GetInverse(cc.normalMass);
        } else {
          cc.pointCount = 1;
        }
      }
    }
  }
  InitVelocityConstraints(step) {
    for (let i = 0; i < this.m_constraintCount; ++i) {
      const c = this.m_constraints[i];
      const bodyA = c.bodyA;
      const bodyB = c.bodyB;
      const invMassA = bodyA.m_invMass;
      const invIA = bodyA.m_invI;
      const invMassB = bodyB.m_invMass;
      const invIB = bodyB.m_invI;
      const normalX = c.normal.x;
      const normalY = c.normal.y;
      const tangentX = normalY;
      const tangentY = -normalX;
      let tCount;
      if (step.warmStarting) {
        tCount = c.pointCount;
        for (let j = 0; j < tCount; ++j) {
          const ccp = c.points[j];
          ccp.normalImpulse *= step.dtRatio;
          ccp.tangentImpulse *= step.dtRatio;
          const PX = ccp.normalImpulse * normalX + ccp.tangentImpulse * tangentX;
          const PY = ccp.normalImpulse * normalY + ccp.tangentImpulse * tangentY;
          bodyA.m_angularVelocity -= invIA * (ccp.rA.x * PY - ccp.rA.y * PX);
          bodyA.m_linearVelocity.x -= invMassA * PX;
          bodyA.m_linearVelocity.y -= invMassA * PY;
          bodyB.m_angularVelocity += invIB * (ccp.rB.x * PY - ccp.rB.y * PX);
          bodyB.m_linearVelocity.x += invMassB * PX;
          bodyB.m_linearVelocity.y += invMassB * PY;
        }
      } else {
        tCount = c.pointCount;
        for (let j = 0; j < tCount; ++j) {
          const ccp2 = c.points[j];
          ccp2.normalImpulse = 0.0;
          ccp2.tangentImpulse = 0.0;
        }
      }
    }
  }
  SolveVelocityConstraints() {
    let ccp;
    let rAX;
    let rAY;
    let rBX;
    let rBY;
    let dvX;
    let dvY;
    let vn;
    let vt;
    let lambda;
    let maxFriction;
    let newImpulse;
    let tangentImpulse;
    let normalImpulse;
    let PX;
    let PY;
    let dX;
    let dY;
    let P1X;
    let P1Y;
    let P2X;
    let P2Y;
    let tMat;
    const constraints = this.m_constraints;
    for (let i = 0; i < this.m_constraintCount; ++i) {
      const c = constraints[i];
      const bodyA = c.bodyA;
      const bodyB = c.bodyB;
      let wA = bodyA.m_angularVelocity;
      let wB = bodyB.m_angularVelocity;
      const vA = bodyA.m_linearVelocity;
      const vB = bodyB.m_linearVelocity;
      let vAX = vA.x;
      let vAY = vA.y;
      let vBX = vB.x;
      let vBY = vB.y;
      const invMassA = bodyA.m_invMass;
      const invIA = bodyA.m_invI;
      const invMassB = bodyB.m_invMass;
      const invIB = bodyB.m_invI;
      const normalX = c.normal.x;
      const normalY = c.normal.y;
      const tangentX = normalY;
      const tangentY = -normalX;
      const friction = c.friction;
      const pointCount = c.pointCount;
      const points = c.points;
      for (let j = 0; j < pointCount; j++) {
        ccp = points[j];
        rAX = ccp.rA.x;
        rAY = ccp.rA.y;
        rBX = ccp.rB.x;
        rBY = ccp.rB.y;
        dvX = vBX - wB * rBY - vAX + wA * rAY;
        dvY = vBY + wB * rBX - vAY - wA * rAX;
        vt = dvX * tangentX + dvY * tangentY;
        lambda = ccp.tangentMass * -vt;
        maxFriction = friction * ccp.normalImpulse;
        tangentImpulse = ccp.tangentImpulse;
        newImpulse = b2Math.Clamp(tangentImpulse + lambda, -maxFriction, maxFriction);
        lambda = newImpulse - tangentImpulse;
        PX = lambda * tangentX;
        PY = lambda * tangentY;
        vAX -= invMassA * PX;
        vAY -= invMassA * PY;
        wA -= invIA * (rAX * PY - rAY * PX);
        vBX += invMassB * PX;
        vBY += invMassB * PY;
        wB += invIB * (rBX * PY - rBY * PX);
        ccp.tangentImpulse = newImpulse;
      }
      if (pointCount == 1) {
        ccp = points[0];
        rAX = ccp.rA.x;
        rAY = ccp.rA.y;
        rBX = ccp.rB.x;
        rBY = ccp.rB.y;
        dvX = vBX + (-wB * rBY) - vAX - (-wA * rAY);
        dvY = vBY + (wB * rBX) - vAY - (wA * rAX);
        vn = dvX * normalX + dvY * normalY;
        normalImpulse = ccp.normalImpulse;
        lambda = -ccp.normalMass * (vn - ccp.velocityBias);
        newImpulse = normalImpulse + lambda;
        newImpulse = newImpulse > 0 ? newImpulse : 0.0;
        lambda = newImpulse - normalImpulse;
        PX = lambda * normalX;
        PY = lambda * normalY;
        vAX -= invMassA * PX;
        vAY -= invMassA * PY;
        wA -= invIA * (rAX * PY - rAY * PX);
        vBX += invMassB * PX;
        vBY += invMassB * PY;
        wB += invIB * (rBX * PY - rBY * PX);
        ccp.normalImpulse = newImpulse;
      } else {
        const cp1 = points[0];
        const cp2 = points[1];
        const rA1X = cp1.rA.x;
        const rA1Y = cp1.rA.y;
        const rB1X = cp1.rB.x;
        const rB1Y = cp1.rB.y;
        const rA2X = cp2.rA.x;
        const rA2Y = cp2.rA.y;
        const rB2X = cp2.rB.x;
        const rB2Y = cp2.rB.y;
        const aX = cp1.normalImpulse;
        const aY = cp2.normalImpulse;
        const dv1X = vBX - wB * rB1Y - vAX + wA * rA1Y;
        const dv1Y = vBY + wB * rB1X - vAY - wA * rA1X;
        const dv2X = vBX - wB * rB2Y - vAX + wA * rA2Y;
        const dv2Y = vBY + wB * rB2X - vAY - wA * rA2X;
        let vn1 = dv1X * normalX + dv1Y * normalY;
        let vn2 = dv2X * normalX + dv2Y * normalY;
        let bX = vn1 - cp1.velocityBias;
        let bY = vn2 - cp2.velocityBias;
        tMat = c.K;
        bX -= tMat.col1.x * aX + tMat.col2.x * aY;
        bY -= tMat.col1.y * aX + tMat.col2.y * aY;
        for (;;) {
          tMat = c.normalMass;
          let xX = -(tMat.col1.x * bX + tMat.col2.x * bY);
          let xY = -(tMat.col1.y * bX + tMat.col2.y * bY);
          if (xX >= 0.0 && xY >= 0.0) {
            dX = xX - aX;
            dY = xY - aY;
            P1X = dX * normalX;
            P1Y = dX * normalY;
            P2X = dY * normalX;
            P2Y = dY * normalY;
            vAX -= invMassA * (P1X + P2X);
            vAY -= invMassA * (P1Y + P2Y);
            wA -= invIA * (rA1X * P1Y - rA1Y * P1X + rA2X * P2Y - rA2Y * P2X);
            vBX += invMassB * (P1X + P2X);
            vBY += invMassB * (P1Y + P2Y);
            wB += invIB * (rB1X * P1Y - rB1Y * P1X + rB2X * P2Y - rB2Y * P2X);
            cp1.normalImpulse = xX;
            cp2.normalImpulse = xY;
            break;
          }
          xX = -cp1.normalMass * bX;
          xY = 0.0;
          vn2 = c.K.col1.y * xX + bY;
          if (xX >= 0.0 && vn2 >= 0.0) {
            dX = xX - aX;
            dY = xY - aY;
            P1X = dX * normalX;
            P1Y = dX * normalY;
            P2X = dY * normalX;
            P2Y = dY * normalY;
            vAX -= invMassA * (P1X + P2X);
            vAY -= invMassA * (P1Y + P2Y);
            wA -= invIA * (rA1X * P1Y - rA1Y * P1X + rA2X * P2Y - rA2Y * P2X);
            vBX += invMassB * (P1X + P2X);
            vBY += invMassB * (P1Y + P2Y);
            wB += invIB * (rB1X * P1Y - rB1Y * P1X + rB2X * P2Y - rB2Y * P2X);
            cp1.normalImpulse = xX;
            cp2.normalImpulse = xY;
            break;
          }
          xX = 0.0;
          xY = -cp2.normalMass * bY;
          vn1 = c.K.col2.x * xY + bX;
          if (xY >= 0.0 && vn1 >= 0.0) {
            dX = xX - aX;
            dY = xY - aY;
            P1X = dX * normalX;
            P1Y = dX * normalY;
            P2X = dY * normalX;
            P2Y = dY * normalY;
            vAX -= invMassA * (P1X + P2X);
            vAY -= invMassA * (P1Y + P2Y);
            wA -= invIA * (rA1X * P1Y - rA1Y * P1X + rA2X * P2Y - rA2Y * P2X);
            vBX += invMassB * (P1X + P2X);
            vBY += invMassB * (P1Y + P2Y);
            wB += invIB * (rB1X * P1Y - rB1Y * P1X + rB2X * P2Y - rB2Y * P2X);
            cp1.normalImpulse = xX;
            cp2.normalImpulse = xY;
            break;
          }
          xX = 0.0;
          xY = 0.0;
          vn1 = bX;
          vn2 = bY;
          if (vn1 >= 0.0 && vn2 >= 0.0) {
            dX = xX - aX;
            dY = xY - aY;
            P1X = dX * normalX;
            P1Y = dX * normalY;
            P2X = dY * normalX;
            P2Y = dY * normalY;
            vAX -= invMassA * (P1X + P2X);
            vAY -= invMassA * (P1Y + P2Y);
            wA -= invIA * (rA1X * P1Y - rA1Y * P1X + rA2X * P2Y - rA2Y * P2X);
            vBX += invMassB * (P1X + P2X);
            vBY += invMassB * (P1Y + P2Y);
            wB += invIB * (rB1X * P1Y - rB1Y * P1X + rB2X * P2Y - rB2Y * P2X);
            cp1.normalImpulse = xX;
            cp2.normalImpulse = xY;
            break;
          }
          break;
        }
      }
      vA.x = vAX;
      vA.y = vAY;
      vB.x = vBX;
      vB.y = vBY;
      bodyA.m_angularVelocity = wA;
      bodyB.m_angularVelocity = wB;
    }
  }
  FinalizeVelocityConstraints() {
    for (let i = 0; i < this.m_constraintCount; ++i) {
      const c = this.m_constraints[i];
      const m = c.manifold;
      for (let j = 0; j < c.pointCount; ++j) {
        const point1 = m.m_points[j];
        const point2 = c.points[j];
        point1.m_normalImpulse = point2.normalImpulse;
        point1.m_tangentImpulse = point2.tangentImpulse;
      }
    }
  }
  SolvePositionConstraints(baumgarte = 0) {
    let minSeparation = 0.0;
    const positionManifold = b2ContactSolver.s_psm;
    const manifoldPoints = positionManifold.m_points;
    const separations = positionManifold.m_separations;
    const constraints = this.m_constraints;
    for (let i = 0; i < this.m_constraintCount; i++) {
      const c = constraints[i];
      const bodyA = c.bodyA;
      const bodyB = c.bodyB;
      const invMassA = bodyA.m_mass * bodyA.m_invMass;
      const invIA = bodyA.m_mass * bodyA.m_invI;
      const invMassB = bodyB.m_mass * bodyB.m_invMass;
      const invIB = bodyB.m_mass * bodyB.m_invI;
      positionManifold.Initialize(c);
      const normalX = positionManifold.m_normal.x;
      const normalY = positionManifold.m_normal.y;
      const sweepA = bodyA.m_sweep;
      const sweepB = bodyB.m_sweep;
      let centerAX = sweepA.c.x;
      let centerAY = sweepA.c.y;
      let angleA = sweepA.a;
      let centerBX = sweepB.c.x;
      let centerBY = sweepB.c.y;
      let angleB = sweepB.a;
      const constraintPoints = c.points;
      const pointCount = c.pointCount;
      for (let j = 0; j < pointCount; j++) {
        const ccp = constraintPoints[j];
        const point = manifoldPoints[j];
        const separation = separations[j];
        const rAX = point.x - centerAX;
        const rAY = point.y - centerAY;
        const rBX = point.x - centerBX;
        const rBY = point.y - centerBY;
        minSeparation = minSeparation < separation ? minSeparation : separation;
        const C = b2Math.Clamp(baumgarte * (separation + b2Settings.b2_linearSlop), -b2Settings.b2_maxLinearCorrection, 0.0);
        const impulse = -ccp.equalizedMass * C;
        const PX = impulse * normalX;
        const PY = impulse * normalY;
        centerAX -= invMassA * PX;
        centerAY -= invMassA * PY;
        angleA -= invIA * (rAX * PY - rAY * PX);
        centerBX += invMassB * PX;
        centerBY += invMassB * PY;
        angleB += invIB * (rBX * PY - rBY * PX);
      }
      if (pointCount > 0) {
        sweepA.c.x = centerAX;
        sweepA.c.y = centerAY;
        sweepA.a = angleA;
        sweepB.c.x = centerBX;
        sweepB.c.y = centerBY;
        sweepB.a = angleB;
        bodyA.SynchronizeTransform();
        bodyB.SynchronizeTransform();
      }
    }
    return minSeparation > (-1.5 * b2Settings.b2_linearSlop);
  }
};

const b2EdgeAndCircleContact = class extends b2Contact {
  static Create(allocator) {
    return new b2EdgeAndCircleContact();
  }
  static Destroy(contact, allocator) {}
  Reset(fixtureA, fixtureB) {
    super.Reset(fixtureA, fixtureB);
  }
  Evaluate() {
    const bA = this.m_fixtureA.GetBody();
    const bB = this.m_fixtureB.GetBody();
    this.b2CollideEdgeAndCircle(this.m_manifold, this.m_fixtureA.GetShape(), bA.m_xf, this.m_fixtureB.GetShape(), bB.m_xf);
  }
  b2CollideEdgeAndCircle(manifold, edge, xf1, circle, xf2) {}
};

const b2NullContact = class extends b2Contact {
  Evaluate() {}
};

const b2PolyAndCircleContact = class extends b2Contact {
  static Create(allocator) {
    return new b2PolyAndCircleContact();
  }
  static Destroy(contact, allocator) {}
  Reset(fixtureA, fixtureB) {
    super.Reset(fixtureA, fixtureB);
    b2Settings.b2Assert(fixtureA.GetType() == b2Shape.e_polygonShape);
    b2Settings.b2Assert(fixtureB.GetType() == b2Shape.e_circleShape);
  }
  Evaluate() {
    const fixtureA = this.m_fixtureA;
    const fixtureB = this.m_fixtureB;
    b2Collision.CollidePolygonAndCircle(this.m_manifold, fixtureA.m_shape, fixtureA.m_body.m_xf, fixtureB.m_shape, fixtureB.m_body.m_xf);
  }
};

const b2PolyAndEdgeContact = class extends b2Contact {
  static Create(allocator) {
    return new b2PolyAndEdgeContact();
  }
  static Destroy(contact, allocator) {}
  Reset(fixtureA, fixtureB) {
    super.Reset(fixtureA, fixtureB);
    b2Settings.b2Assert(fixtureA.GetType() == b2Shape.e_polygonShape);
    b2Settings.b2Assert(fixtureB.GetType() == b2Shape.e_edgeShape);
  }
  Evaluate() {
    const bA = this.m_fixtureA.GetBody();
    const bB = this.m_fixtureB.GetBody();
    this.b2CollidePolyAndEdge(this.m_manifold, this.m_fixtureA.GetShape(), bA.m_xf, this.m_fixtureB.GetShape(), bB.m_xf);
  }
  b2CollidePolyAndEdge(manifold, polygon, xf1, edge, xf2) {}
};

const b2PolygonContact = class extends b2Contact {
  static Create(allocator) {
    return new b2PolygonContact();
  }
  static Destroy(contact, allocator) {}
  Reset(fixtureA, fixtureB) {
    super.Reset(fixtureA, fixtureB);
  }
  Evaluate() {
    const fixtureA = this.m_fixtureA;
    const fixtureB = this.m_fixtureB;
    b2Collision.CollidePolygons(this.m_manifold, fixtureA.m_shape, fixtureA.m_body.m_xf, fixtureB.m_shape, fixtureB.m_body.m_xf);
  }
};

const b2Controller = class {
  Step(step) {}
  Draw(debugDraw) {}
  AddBody(body) {
    const edge = new b2ControllerEdge();
    edge.controller = this;
    edge.body = body;
    edge.nextBody = this.m_bodyList;
    edge.prevBody = null;
    this.m_bodyList = edge;
    if (edge.nextBody) edge.nextBody.prevBody = edge;
    this.m_bodyCount++;
    edge.nextController = body.m_controllerList;
    edge.prevController = null;
    body.m_controllerList = edge;
    if (edge.nextController) edge.nextController.prevController = edge;
    body.m_controllerCount++;
  }
  RemoveBody(body) {
    let edge = body.m_controllerList;
    while (edge && edge.controller != this)
    edge = edge.nextController;
    if (edge.prevBody) edge.prevBody.nextBody = edge.nextBody;
    if (edge.nextBody) edge.nextBody.prevBody = edge.prevBody;
    if (edge.nextController) edge.nextController.prevController = edge.prevController;
    if (edge.prevController) edge.prevController.nextController = edge.nextController;
    if (this.m_bodyList == edge) this.m_bodyList = edge.nextBody;
    if (body.m_controllerList == edge) body.m_controllerList = edge.nextController;
    body.m_controllerCount--;
    this.m_bodyCount--;
  }
  Clear() {
    while (this.m_bodyList)
    this.RemoveBody(this.m_bodyList.body);
  }
  GetNext() {
    return this.m_next;
  }
  GetWorld() {
    return this.m_world;
  }
  GetBodyList() {
    return this.m_bodyList;
  }
};

const b2BuoyancyController = class extends b2Controller {
  constructor() {
    super();
    this.normal = new b2Vec2(0, -1);
    this.offset = 0;
    this.density = 0;
    this.velocity = new b2Vec2(0, 0);
    this.linearDrag = 2;
    this.angularDrag = 1;
    this.useDensity = false;
    this.useWorldGravity = true;
    this.gravity = null;
  }
  Step(step) {
    if (!this.m_bodyList) return;
    if (this.useWorldGravity) {
      this.gravity = this.GetWorld().GetGravity().Copy();
    }
    for (let i = this.m_bodyList; i; i = i.nextBody) {
      const body = i.body;
      if (body.IsAwake() == false) {
        continue;
      }
      const areac = new b2Vec2();
      const massc = new b2Vec2();
      let area = 0.0;
      let mass = 0.0;
      for (let fixture = body.GetFixtureList(); fixture; fixture = fixture.GetNext()) {
        const sc = new b2Vec2();
        const sarea = fixture.GetShape().ComputeSubmergedArea(this.normal, this.offset, body.GetTransform(), sc);
        area += sarea;
        areac.x += sarea * sc.x;
        areac.y += sarea * sc.y;
        let shapeDensity;
        if (this.useDensity) {
          shapeDensity = 1;
        } else {
          shapeDensity = 1;
        }
        mass += sarea * shapeDensity;
        massc.x += sarea * sc.x * shapeDensity;
        massc.y += sarea * sc.y * shapeDensity;
      }
      areac.x /= area;
      areac.y /= area;
      massc.x /= mass;
      massc.y /= mass;
      if (area < Number.MIN_VALUE) continue;
      const buoyancyForce = this.gravity.GetNegative();
      buoyancyForce.Multiply(this.density * area);
      body.ApplyForce(buoyancyForce, massc);
      const dragForce = body.GetLinearVelocityFromWorldPoint(areac);
      dragForce.Subtract(this.velocity);
      dragForce.Multiply(-this.linearDrag * area);
      body.ApplyForce(dragForce, areac);
      body.ApplyTorque(-body.GetInertia() / body.GetMass() * area * body.GetAngularVelocity() * this.angularDrag);
    }
  }
  Draw(debugDraw) {
    const r = 1000;
    const p1 = new b2Vec2();
    const p2 = new b2Vec2();
    p1.x = this.normal.x * this.offset + this.normal.y * r;
    p1.y = this.normal.y * this.offset - this.normal.x * r;
    p2.x = this.normal.x * this.offset - this.normal.y * r;
    p2.y = this.normal.y * this.offset + this.normal.x * r;
    const color = new b2Color(0, 0, 1);
    debugDraw.DrawSegment(p1, p2, color);
  }
};

const b2ConstantAccelController = class extends b2Controller {
  constructor() {
    super();
    this.A = new b2Vec2(0, 0);
  }
  Step(step) {
    const smallA = new b2Vec2(this.A.x * step.dt, this.A.y * step.dt);
    for (let i = this.m_bodyList; i; i = i.nextBody) {
      const body = i.body;
      if (!body.IsAwake()) continue;
      body.SetLinearVelocity(new b2Vec2(body.GetLinearVelocity().x + smallA.x, body.GetLinearVelocity().y + smallA.y));
    }
  }
};

const b2ConstantForceController = class extends b2Controller {
  constructor() {
    super();
    this.F = new b2Vec2(0, 0);
  }
  Step(step) {
    for (let i = this.m_bodyList; i; i = i.nextBody) {
      const body = i.body;
      if (!body.IsAwake()) continue;
      body.ApplyForce(this.F, body.GetWorldCenter());
    }
  }
};

const b2ControllerEdge = class {
};

const b2GravityController = class extends b2Controller {
  constructor() {
    super();
    this.G = 1;
    this.invSqr = true;
  }
  Step(step) {
    let body1;
    let p1;
    let mass1;
    let body2;
    let p2;
    let dx;
    let dy;
    let r2;
    let f;
    if (this.invSqr) {
      for (let i = this.m_bodyList; i; i = i.nextBody) {
        body1 = i.body;
        p1 = body1.GetWorldCenter();
        mass1 = body1.GetMass();
        for (let j = this.m_bodyList; j != i; j = j.nextBody) {
          body2 = j.body;
          p2 = body2.GetWorldCenter();
          dx = p2.x - p1.x;
          dy = p2.y - p1.y;
          r2 = dx * dx + dy * dy;
          if (r2 < Number.MIN_VALUE) continue;
          f = new b2Vec2(dx, dy);
          f.Multiply(this.G / r2 / Math.sqrt(r2) * mass1 * body2.GetMass());
          if (body1.IsAwake()) body1.ApplyForce(f, p1);
          f.Multiply(-1);
          if (body2.IsAwake()) body2.ApplyForce(f, p2);
        }
      }
    } else {
      for (let i = this.m_bodyList; i; i = i.nextBody) {
        body1 = i.body;
        p1 = body1.GetWorldCenter();
        mass1 = body1.GetMass();
        for (let j = this.m_bodyList; j != i; j = j.nextBody) {
          body2 = j.body;
          p2 = body2.GetWorldCenter();
          dx = p2.x - p1.x;
          dy = p2.y - p1.y;
          r2 = dx * dx + dy * dy;
          if (r2 < Number.MIN_VALUE) continue;
          f = new b2Vec2(dx, dy);
          f.Multiply(this.G / r2 * mass1 * body2.GetMass());
          if (body1.IsAwake()) body1.ApplyForce(f, p1);
          f.Multiply(-1);
          if (body2.IsAwake()) body2.ApplyForce(f, p2);
        }
      }
    }
  }
};

const b2TensorDampingController = class extends b2Controller {
  constructor() {
    super();
    this.T = new b2Mat22();
    this.maxTimestep = 0;
  }
  SetAxisAligned(xDamping = 0, yDamping = 0) {
    this.T.col1.x = -xDamping;
    this.T.col1.y = 0;
    this.T.col2.x = 0;
    this.T.col2.y = -yDamping;
    if (xDamping > 0 || yDamping > 0) {
      this.maxTimestep = 1 / Math.max(xDamping, yDamping);
    } else {
      this.maxTimestep = 0;
    }
  }
  Step(step) {
    let timestep = step.dt;
    if (timestep <= Number.MIN_VALUE) return;
    if (timestep > this.maxTimestep && this.maxTimestep > 0) timestep = this.maxTimestep;
    for (let i = this.m_bodyList; i; i = i.nextBody) {
      const body = i.body;
      if (!body.IsAwake()) {
        continue;
      }
      const damping = body.GetWorldVector(b2Math.MulMV(this.T, body.GetLocalVector(body.GetLinearVelocity())));
      body.SetLinearVelocity(new b2Vec2(body.GetLinearVelocity().x + damping.x * timestep, body.GetLinearVelocity().y + damping.y * timestep));
    }
  }
};

const b2Joint = class {
  static e_unknownJoint = 0;
  static e_revoluteJoint = 1;
  static e_prismaticJoint = 2;
  static e_distanceJoint = 3;
  static e_pulleyJoint = 4;
  static e_mouseJoint = 5;
  static e_gearJoint = 6;
  static e_lineJoint = 7;
  static e_weldJoint = 8;
  static e_frictionJoint = 9;
  static e_inactiveLimit = 0;
  static e_atLowerLimit = 1;
  static e_atUpperLimit = 2;
  static e_equalLimits = 3;
  constructor(def) {
    this.m_edgeA = new b2JointEdge();
    this.m_edgeB = new b2JointEdge();
    this.m_localCenterA = new b2Vec2();
    this.m_localCenterB = new b2Vec2();
    b2Settings.b2Assert(def.bodyA != def.bodyB);
    this.m_type = def.type;
    this.m_prev = null;
    this.m_next = null;
    this.m_bodyA = def.bodyA;
    this.m_bodyB = def.bodyB;
    this.m_collideConnected = def.collideConnected;
    this.m_islandFlag = false;
    this.m_userData = def.userData;
  }
  GetType() {
    return this.m_type;
  }
  GetAnchorA() {
    return null;
  }
  GetAnchorB() {
    return null;
  }
  GetReactionForce(inv_dt) {
    return null;
  }
  GetReactionTorque(inv_dt) {
    return 0.0;
  }
  GetBodyA() {
    return this.m_bodyA;
  }
  GetBodyB() {
    return this.m_bodyB;
  }
  GetNext() {
    return this.m_next;
  }
  GetUserData() {
    return this.m_userData;
  }
  SetUserData(data) {
    this.m_userData = data;
  }
  IsActive() {
    return this.m_bodyA.IsActive() && this.m_bodyB.IsActive();
  }
  static Create(def, allocator) {
    let joint = null;
    switch (def.type) {
    case b2Joint.e_distanceJoint:
      joint = new b2DistanceJoint(def);
      break;
    case b2Joint.e_mouseJoint:
      joint = new b2MouseJoint(def);
      break;
    case b2Joint.e_prismaticJoint:
      joint = new b2PrismaticJoint(def);
      break;
    case b2Joint.e_revoluteJoint:
      joint = new b2RevoluteJoint(def);
      break;
    case b2Joint.e_pulleyJoint:
      joint = new b2PulleyJoint(def);
      break;
    case b2Joint.e_gearJoint:
      joint = new b2GearJoint(def);
      break;
    case b2Joint.e_lineJoint:
      joint = new b2LineJoint(def);
      break;
    case b2Joint.e_weldJoint:
      joint = new b2WeldJoint(def);
      break;
    case b2Joint.e_frictionJoint:
      joint = new b2FrictionJoint(def);
      break;
    default:
      break;
    }
    return joint;
  }
  static Destroy(joint, allocator) {}
  InitVelocityConstraints(step) {}
  SolveVelocityConstraints(step) {}
  FinalizeVelocityConstraints() {}
  SolvePositionConstraints(baumgarte) {
    return false;
  }
};

const b2DistanceJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_u = new b2Vec2();
    this.m_localAnchor1.SetV(def.localAnchorA);
    this.m_localAnchor2.SetV(def.localAnchorB);
    this.m_length = def.length;
    this.m_frequencyHz = def.frequencyHz;
    this.m_dampingRatio = def.dampingRatio;
    this.m_impulse = 0.0;
    this.m_gamma = 0.0;
    this.m_bias = 0.0;
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse * this.m_u.x, inv_dt * this.m_impulse * this.m_u.y);
  }
  GetReactionTorque(inv_dt) {
    return 0.0;
  }
  GetLength() {
    return this.m_length;
  }
  SetLength(length = 0) {
    this.m_length = length;
  }
  GetFrequency() {
    return this.m_frequencyHz;
  }
  SetFrequency(hz = 0) {
    this.m_frequencyHz = hz;
  }
  GetDampingRatio() {
    return this.m_dampingRatio;
  }
  SetDampingRatio(ratio = 0) {
    this.m_dampingRatio = ratio;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    this.m_u.x = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
    this.m_u.y = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
    const length = Math.sqrt(this.m_u.x * this.m_u.x + this.m_u.y * this.m_u.y);
    if (length > b2Settings.b2_linearSlop) {
      this.m_u.Multiply(1.0 / length);
    } else {
      this.m_u.SetZero();
    }
    const cr1u = r1X * this.m_u.y - r1Y * this.m_u.x;
    const cr2u = r2X * this.m_u.y - r2Y * this.m_u.x;
    const invMass = bA.m_invMass + bA.m_invI * cr1u * cr1u + bB.m_invMass + bB.m_invI * cr2u * cr2u;
    this.m_mass = invMass != 0.0 ? 1.0 / invMass : 0.0;
    if (this.m_frequencyHz > 0.0) {
      const C = length - this.m_length;
      const omega = 2.0 * Math.PI * this.m_frequencyHz;
      const d = 2.0 * this.m_mass * this.m_dampingRatio * omega;
      const k = this.m_mass * omega * omega;
      this.m_gamma = step.dt * (d + step.dt * k);
      this.m_gamma = this.m_gamma != 0.0 ? 1 / this.m_gamma : 0.0;
      this.m_bias = C * step.dt * k * this.m_gamma;
      this.m_mass = invMass + this.m_gamma;
      this.m_mass = this.m_mass != 0.0 ? 1.0 / this.m_mass : 0.0;
    }
    if (step.warmStarting) {
      this.m_impulse *= step.dtRatio;
      const PX = this.m_impulse * this.m_u.x;
      const PY = this.m_impulse * this.m_u.y;
      bA.m_linearVelocity.x -= bA.m_invMass * PX;
      bA.m_linearVelocity.y -= bA.m_invMass * PY;
      bA.m_angularVelocity -= bA.m_invI * (r1X * PY - r1Y * PX);
      bB.m_linearVelocity.x += bB.m_invMass * PX;
      bB.m_linearVelocity.y += bB.m_invMass * PY;
      bB.m_angularVelocity += bB.m_invI * (r2X * PY - r2Y * PX);
    } else {
      this.m_impulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const v1X = bA.m_linearVelocity.x + (-bA.m_angularVelocity * r1Y);
    const v1Y = bA.m_linearVelocity.y + (bA.m_angularVelocity * r1X);
    const v2X = bB.m_linearVelocity.x + (-bB.m_angularVelocity * r2Y);
    const v2Y = bB.m_linearVelocity.y + (bB.m_angularVelocity * r2X);
    const Cdot = this.m_u.x * (v2X - v1X) + this.m_u.y * (v2Y - v1Y);
    const impulse = -this.m_mass * (Cdot + this.m_bias + this.m_gamma * this.m_impulse);
    this.m_impulse += impulse;
    const PX = impulse * this.m_u.x;
    const PY = impulse * this.m_u.y;
    bA.m_linearVelocity.x -= bA.m_invMass * PX;
    bA.m_linearVelocity.y -= bA.m_invMass * PY;
    bA.m_angularVelocity -= bA.m_invI * (r1X * PY - r1Y * PX);
    bB.m_linearVelocity.x += bB.m_invMass * PX;
    bB.m_linearVelocity.y += bB.m_invMass * PY;
    bB.m_angularVelocity += bB.m_invI * (r2X * PY - r2Y * PX);
  }
  SolvePositionConstraints(baumgarte) {
    if (this.m_frequencyHz > 0.0) {
      return true;
    }
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    let dX = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
    let dY = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
    let length = Math.sqrt(dX * dX + dY * dY);
    if (length < b2Settings.b2_epsilon) {
      length = 0.0;
    } else {
      const invLength = 1.0 / length;
      dX *= invLength;
      dY *= invLength;
    }
    let C = length - this.m_length;
    C = b2Math.Clamp(C, -b2Settings.b2_maxLinearCorrection, b2Settings.b2_maxLinearCorrection);
    const impulse = -this.m_mass * C;
    this.m_u.Set(dX, dY);
    const PX = impulse * this.m_u.x;
    const PY = impulse * this.m_u.y;
    bA.m_sweep.c.x -= bA.m_invMass * PX;
    bA.m_sweep.c.y -= bA.m_invMass * PY;
    bA.m_sweep.a -= bA.m_invI * (r1X * PY - r1Y * PX);
    bB.m_sweep.c.x += bB.m_invMass * PX;
    bB.m_sweep.c.y += bB.m_invMass * PY;
    bB.m_sweep.a += bB.m_invI * (r2X * PY - r2Y * PX);
    bA.SynchronizeTransform();
    bB.SynchronizeTransform();
    return b2Math.Abs(C) < b2Settings.b2_linearSlop;
  }
};

const b2JointDef = class {
  constructor() {
    this.type = b2Joint.e_unknownJoint;
    this.userData = null;
    this.bodyA = null;
    this.bodyB = null;
    this.collideConnected = false;
  }
};

const b2DistanceJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.type = b2Joint.e_distanceJoint;
    this.length = 1.0;
    this.frequencyHz = 0.0;
    this.dampingRatio = 0.0;
  }
  Initialize(bA, bB, anchorA, anchorB) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA.SetV(this.bodyA.GetLocalPoint(anchorA));
    this.localAnchorB.SetV(this.bodyB.GetLocalPoint(anchorB));
    const dX = anchorB.x - anchorA.x;
    const dY = anchorB.y - anchorA.y;
    this.length = Math.sqrt(dX * dX + dY * dY);
    this.frequencyHz = 0.0;
    this.dampingRatio = 0.0;
  }
};

const b2FrictionJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_localAnchorA = new b2Vec2();
    this.m_localAnchorB = new b2Vec2();
    this.m_linearMass = new b2Mat22();
    this.m_linearImpulse = new b2Vec2();
    this.m_tempK = new b2Mat22();
    this.m_localAnchorA.SetV(def.localAnchorA);
    this.m_localAnchorB.SetV(def.localAnchorB);
    this.m_linearMass.SetZero();
    this.m_angularMass = 0.0;
    this.m_linearImpulse.SetZero();
    this.m_angularImpulse = 0.0;
    this.m_maxForce = def.maxForce;
    this.m_maxTorque = def.maxTorque;
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchorA);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchorB);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_linearImpulse.x, inv_dt * this.m_linearImpulse.y);
  }
  GetReactionTorque(inv_dt = 0) {
    return inv_dt * this.m_angularImpulse;
  }
  SetMaxForce(force = 0) {
    this.m_maxForce = force;
  }
  GetMaxForce() {
    return this.m_maxForce;
  }
  SetMaxTorque(torque = 0) {
    this.m_maxTorque = torque;
  }
  GetMaxTorque() {
    return this.m_maxTorque;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let rAX = this.m_localAnchorA.x - bA.m_sweep.localCenter.x;
    let rAY = this.m_localAnchorA.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rAX + tMat.col2.x * rAY;
    rAY = tMat.col1.y * rAX + tMat.col2.y * rAY;
    rAX = tX;
    tMat = bB.m_xf.R;
    let rBX = this.m_localAnchorB.x - bB.m_sweep.localCenter.x;
    let rBY = this.m_localAnchorB.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * rBX + tMat.col2.x * rBY;
    rBY = tMat.col1.y * rBX + tMat.col2.y * rBY;
    rBX = tX;
    const mA = bA.m_invMass;
    const mB = bB.m_invMass;
    const iA = bA.m_invI;
    const iB = bB.m_invI;
    const K = this.m_tempK;
    K.col1.x = mA + mB;
    K.col2.x = 0.0;
    K.col1.y = 0.0;
    K.col2.y = mA + mB;
    K.col1.x += iA * rAY * rAY;
    K.col2.x += (-iA * rAX * rAY);
    K.col1.y += (-iA * rAX * rAY);
    K.col2.y += iA * rAX * rAX;
    K.col1.x += iB * rBY * rBY;
    K.col2.x += (-iB * rBX * rBY);
    K.col1.y += (-iB * rBX * rBY);
    K.col2.y += iB * rBX * rBX;
    K.GetInverse(this.m_linearMass);
    this.m_angularMass = iA + iB;
    if (this.m_angularMass > 0.0) {
      this.m_angularMass = 1.0 / this.m_angularMass;
    }
    if (step.warmStarting) {
      this.m_linearImpulse.x *= step.dtRatio;
      this.m_linearImpulse.y *= step.dtRatio;
      this.m_angularImpulse *= step.dtRatio;
      const P = this.m_linearImpulse;
      bA.m_linearVelocity.x -= mA * P.x;
      bA.m_linearVelocity.y -= mA * P.y;
      bA.m_angularVelocity -= iA * (rAX * P.y - rAY * P.x + this.m_angularImpulse);
      bB.m_linearVelocity.x += mB * P.x;
      bB.m_linearVelocity.y += mB * P.y;
      bB.m_angularVelocity += iB * (rBX * P.y - rBY * P.x + this.m_angularImpulse);
    } else {
      this.m_linearImpulse.SetZero();
      this.m_angularImpulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const vA = bA.m_linearVelocity;
    let wA = bA.m_angularVelocity;
    const vB = bB.m_linearVelocity;
    let wB = bB.m_angularVelocity;
    const mA = bA.m_invMass;
    const mB = bB.m_invMass;
    const iA = bA.m_invI;
    const iB = bB.m_invI;
    let tMat = bA.m_xf.R;
    let rAX = this.m_localAnchorA.x - bA.m_sweep.localCenter.x;
    let rAY = this.m_localAnchorA.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rAX + tMat.col2.x * rAY;
    rAY = tMat.col1.y * rAX + tMat.col2.y * rAY;
    rAX = tX;
    tMat = bB.m_xf.R;
    let rBX = this.m_localAnchorB.x - bB.m_sweep.localCenter.x;
    let rBY = this.m_localAnchorB.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * rBX + tMat.col2.x * rBY;
    rBY = tMat.col1.y * rBX + tMat.col2.y * rBY;
    rBX = tX;
    let maxImpulse;
    {
      const Cdot = wB - wA;
      let impulse = -this.m_angularMass * Cdot;
      const oldImpulse = this.m_angularImpulse;
      maxImpulse = step.dt * this.m_maxTorque;
      this.m_angularImpulse = b2Math.Clamp(this.m_angularImpulse + impulse, -maxImpulse, maxImpulse);
      impulse = this.m_angularImpulse - oldImpulse;
      wA -= iA * impulse;
      wB += iB * impulse;
    }
    {
      const CdotX = vB.x - wB * rBY - vA.x + wA * rAY;
      const CdotY = vB.y + wB * rBX - vA.y - wA * rAX;
      tMat = this.m_linearMass;
      let impulseVX = tMat.col1.x * -CdotX + tMat.col2.x * -CdotY;
      let impulseVY = tMat.col1.y * -CdotX + tMat.col2.y * -CdotY;
      const oldImpulseVX = this.m_linearImpulse.x;
      const oldImpulseVY = this.m_linearImpulse.y;
      this.m_linearImpulse.x += impulseVX;
      this.m_linearImpulse.y += impulseVY;
      maxImpulse = step.dt * this.m_maxForce;
      if (this.m_linearImpulse.LengthSquared() > maxImpulse * maxImpulse) {
        this.m_linearImpulse.Normalize();
        this.m_linearImpulse.Multiply(maxImpulse);
      }
      impulseVX = this.m_linearImpulse.x - oldImpulseVX;
      impulseVY = this.m_linearImpulse.y - oldImpulseVY;
      vA.x -= mA * impulseVX;
      vA.y -= mA * impulseVY;
      wA -= iA * (rAX * impulseVY - rAY * impulseVX);
      vB.x += mB * impulseVX;
      vB.y += mB * impulseVY;
      wB += iB * (rBX * impulseVY - rBY * impulseVX);
    }
    bA.m_angularVelocity = wA;
    bB.m_angularVelocity = wB;
  }
  SolvePositionConstraints(baumgarte) {
    return true;
  }
};

const b2FrictionJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.type = b2Joint.e_frictionJoint;
    this.maxForce = 0.0;
    this.maxTorque = 0.0;
  }
  Initialize(bA, bB, anchor) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA.SetV(this.bodyA.GetLocalPoint(anchor));
    this.localAnchorB.SetV(this.bodyB.GetLocalPoint(anchor));
  }
};

const b2GearJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_groundAnchor1 = new b2Vec2();
    this.m_groundAnchor2 = new b2Vec2();
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_J = new b2Jacobian();
    const type1 = def.joint1.m_type;
    const type2 = def.joint2.m_type;
    this.m_revolute1 = null;
    this.m_prismatic1 = null;
    this.m_revolute2 = null;
    this.m_prismatic2 = null;
    let coordinate1;
    let coordinate2;
    this.m_ground1 = def.joint1.GetBodyA();
    this.m_bodyA = def.joint1.GetBodyB();
    if (type1 == b2Joint.e_revoluteJoint) {
      this.m_revolute1 = def.joint1;
      this.m_groundAnchor1.SetV(this.m_revolute1.m_localAnchor1);
      this.m_localAnchor1.SetV(this.m_revolute1.m_localAnchor2);
      coordinate1 = this.m_revolute1.GetJointAngle();
    } else {
      this.m_prismatic1 = def.joint1;
      this.m_groundAnchor1.SetV(this.m_prismatic1.m_localAnchor1);
      this.m_localAnchor1.SetV(this.m_prismatic1.m_localAnchor2);
      coordinate1 = this.m_prismatic1.GetJointTranslation();
    }
    this.m_ground2 = def.joint2.GetBodyA();
    this.m_bodyB = def.joint2.GetBodyB();
    if (type2 == b2Joint.e_revoluteJoint) {
      this.m_revolute2 = def.joint2;
      this.m_groundAnchor2.SetV(this.m_revolute2.m_localAnchor1);
      this.m_localAnchor2.SetV(this.m_revolute2.m_localAnchor2);
      coordinate2 = this.m_revolute2.GetJointAngle();
    } else {
      this.m_prismatic2 = def.joint2;
      this.m_groundAnchor2.SetV(this.m_prismatic2.m_localAnchor1);
      this.m_localAnchor2.SetV(this.m_prismatic2.m_localAnchor2);
      coordinate2 = this.m_prismatic2.GetJointTranslation();
    }
    this.m_ratio = def.ratio;
    this.m_constant = coordinate1 + this.m_ratio * coordinate2;
    this.m_impulse = 0.0;
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse * this.m_J.linearB.x, inv_dt * this.m_impulse * this.m_J.linearB.y);
  }
  GetReactionTorque(inv_dt = 0) {
    const tMat = this.m_bodyB.m_xf.R;
    let rX = this.m_localAnchor2.x - this.m_bodyB.m_sweep.localCenter.x;
    let rY = this.m_localAnchor2.y - this.m_bodyB.m_sweep.localCenter.y;
    const tX = tMat.col1.x * rX + tMat.col2.x * rY;
    rY = tMat.col1.y * rX + tMat.col2.y * rY;
    rX = tX;
    const PX = this.m_impulse * this.m_J.linearB.x;
    const PY = this.m_impulse * this.m_J.linearB.y;
    return inv_dt * (this.m_impulse * this.m_J.angularB - rX * PY + rY * PX);
  }
  GetRatio() {
    return this.m_ratio;
  }
  SetRatio(ratio = 0) {
    this.m_ratio = ratio;
  }
  InitVelocityConstraints(step) {
    const g1 = this.m_ground1;
    const g2 = this.m_ground2;
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let ugX;
    let ugY;
    let rX;
    let rY;
    let tMat;
    let tVec;
    let crug;
    let tX;
    let K = 0.0;
    this.m_J.SetZero();
    if (this.m_revolute1) {
      this.m_J.angularA = -1.0;
      K += bA.m_invI;
    } else {
      tMat = g1.m_xf.R;
      tVec = this.m_prismatic1.m_localXAxis1;
      ugX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      ugY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tMat = bA.m_xf.R;
      rX = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      rY = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * rX + tMat.col2.x * rY;
      rY = tMat.col1.y * rX + tMat.col2.y * rY;
      rX = tX;
      crug = rX * ugY - rY * ugX;
      this.m_J.linearA.Set(-ugX, -ugY);
      this.m_J.angularA = -crug;
      K += bA.m_invMass + bA.m_invI * crug * crug;
    }
    if (this.m_revolute2) {
      this.m_J.angularB = -this.m_ratio;
      K += this.m_ratio * this.m_ratio * bB.m_invI;
    } else {
      tMat = g2.m_xf.R;
      tVec = this.m_prismatic2.m_localXAxis1;
      ugX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      ugY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      tMat = bB.m_xf.R;
      rX = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      rY = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * rX + tMat.col2.x * rY;
      rY = tMat.col1.y * rX + tMat.col2.y * rY;
      rX = tX;
      crug = rX * ugY - rY * ugX;
      this.m_J.linearB.Set(-this.m_ratio * ugX, -this.m_ratio * ugY);
      this.m_J.angularB = -this.m_ratio * crug;
      K += this.m_ratio * this.m_ratio * (bB.m_invMass + bB.m_invI * crug * crug);
    }
    this.m_mass = K > 0.0 ? 1.0 / K : 0.0;
    if (step.warmStarting) {
      bA.m_linearVelocity.x += bA.m_invMass * this.m_impulse * this.m_J.linearA.x;
      bA.m_linearVelocity.y += bA.m_invMass * this.m_impulse * this.m_J.linearA.y;
      bA.m_angularVelocity += bA.m_invI * this.m_impulse * this.m_J.angularA;
      bB.m_linearVelocity.x += bB.m_invMass * this.m_impulse * this.m_J.linearB.x;
      bB.m_linearVelocity.y += bB.m_invMass * this.m_impulse * this.m_J.linearB.y;
      bB.m_angularVelocity += bB.m_invI * this.m_impulse * this.m_J.angularB;
    } else {
      this.m_impulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const Cdot = this.m_J.Compute(bA.m_linearVelocity, bA.m_angularVelocity, bB.m_linearVelocity, bB.m_angularVelocity);
    const impulse = -this.m_mass * Cdot;
    this.m_impulse += impulse;
    bA.m_linearVelocity.x += bA.m_invMass * impulse * this.m_J.linearA.x;
    bA.m_linearVelocity.y += bA.m_invMass * impulse * this.m_J.linearA.y;
    bA.m_angularVelocity += bA.m_invI * impulse * this.m_J.angularA;
    bB.m_linearVelocity.x += bB.m_invMass * impulse * this.m_J.linearB.x;
    bB.m_linearVelocity.y += bB.m_invMass * impulse * this.m_J.linearB.y;
    bB.m_angularVelocity += bB.m_invI * impulse * this.m_J.angularB;
  }
  SolvePositionConstraints(baumgarte) {
    const linearError = 0.0;
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let coordinate1;
    let coordinate2;
    if (this.m_revolute1) {
      coordinate1 = this.m_revolute1.GetJointAngle();
    } else {
      coordinate1 = this.m_prismatic1.GetJointTranslation();
    }
    if (this.m_revolute2) {
      coordinate2 = this.m_revolute2.GetJointAngle();
    } else {
      coordinate2 = this.m_prismatic2.GetJointTranslation();
    }
    const C = this.m_constant - (coordinate1 + this.m_ratio * coordinate2);
    const impulse = -this.m_mass * C;
    bA.m_sweep.c.x += bA.m_invMass * impulse * this.m_J.linearA.x;
    bA.m_sweep.c.y += bA.m_invMass * impulse * this.m_J.linearA.y;
    bA.m_sweep.a += bA.m_invI * impulse * this.m_J.angularA;
    bB.m_sweep.c.x += bB.m_invMass * impulse * this.m_J.linearB.x;
    bB.m_sweep.c.y += bB.m_invMass * impulse * this.m_J.linearB.y;
    bB.m_sweep.a += bB.m_invI * impulse * this.m_J.angularB;
    bA.SynchronizeTransform();
    bB.SynchronizeTransform();
    return linearError < b2Settings.b2_linearSlop;
  }
};

const b2GearJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.type = b2Joint.e_gearJoint;
    this.joint1 = null;
    this.joint2 = null;
    this.ratio = 1.0;
  }
};

const b2Jacobian = class {
  constructor() {
    this.linearA = new b2Vec2();
    this.linearB = new b2Vec2();
  }
  SetZero() {
    this.linearA.SetZero();
    this.angularA = 0.0;
    this.linearB.SetZero();
    this.angularB = 0.0;
  }
  Set(x1, a1 = 0, x2, a2 = 0) {
    this.linearA.SetV(x1);
    this.angularA = a1;
    this.linearB.SetV(x2);
    this.angularB = a2;
  }
  Compute(x1, a1 = 0, x2, a2 = 0) {
    return (this.linearA.x * x1.x + this.linearA.y * x1.y) + this.angularA * a1 + (this.linearB.x * x2.x + this.linearB.y * x2.y) + this.angularB * a2;
  }
};

const b2JointEdge = class {
};

const b2LineJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_localXAxis1 = new b2Vec2();
    this.m_localYAxis1 = new b2Vec2();
    this.m_axis = new b2Vec2();
    this.m_perp = new b2Vec2();
    this.m_K = new b2Mat22();
    this.m_impulse = new b2Vec2();
    this.m_tempImpulse = new b2Vec2();
    this.m_tempRotation1 = new b2Mat22();
    this.m_tempRotation2 = new b2Mat22();
    this.m_localAnchor1.SetV(def.localAnchorA);
    this.m_localAnchor2.SetV(def.localAnchorB);
    this.m_localXAxis1.SetV(def.localAxisA);
    this.m_localYAxis1.x = -this.m_localXAxis1.y;
    this.m_localYAxis1.y = this.m_localXAxis1.x;
    this.m_impulse.SetZero();
    this.m_motorMass = 0.0;
    this.m_motorImpulse = 0.0;
    this.m_lowerTranslation = def.lowerTranslation;
    this.m_upperTranslation = def.upperTranslation;
    this.m_maxMotorForce = def.maxMotorForce;
    this.m_motorSpeed = def.motorSpeed;
    this.m_enableLimit = def.enableLimit;
    this.m_enableMotor = def.enableMotor;
    this.m_limitState = b2Joint.e_inactiveLimit;
    this.m_axis.SetZero();
    this.m_perp.SetZero();
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * (this.m_impulse.x * this.m_perp.x + (this.m_motorImpulse + this.m_impulse.y) * this.m_axis.x), inv_dt * (this.m_impulse.x * this.m_perp.y + (this.m_motorImpulse + this.m_impulse.y) * this.m_axis.y));
  }
  GetReactionTorque(inv_dt) {
    return 0.0;
  }
  GetJointTranslation() {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let tVec = this.m_localAnchor1;
    const p1X = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + bA.m_xf.position.x;
    const p1Y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + bA.m_xf.position.y;
    tVec = this.m_localXAxis1;
    const axisX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    const axisY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    tMat = bB.m_xf.R;
    tVec = this.m_localAnchor2;
    const p2X = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + bB.m_xf.position.x;
    const p2Y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + bB.m_xf.position.y;
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    const translation = axisX * dX + axisY * dY;
    return translation;
  }
  GetJointSpeed() {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const p1X = bA.m_sweep.c.x + r1X;
    const p1Y = bA.m_sweep.c.y + r1Y;
    const p2X = bB.m_sweep.c.x + r2X;
    const p2Y = bB.m_sweep.c.y + r2Y;
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    const axis = bA.GetWorldVector(this.m_localXAxis1);
    const v1 = bA.m_linearVelocity;
    const v2 = bB.m_linearVelocity;
    const w1 = bA.m_angularVelocity;
    const w2 = bB.m_angularVelocity;
    const speed = (dX * (-w1 * axis.y) + dY * (w1 * axis.x)) + (axis.x * (((v2.x + (-w2 * r2Y)) - v1.x) - (-w1 * r1Y)) + axis.y * (((v2.y + (w2 * r2X)) - v1.y) - (w1 * r1X)));
    return speed;
  }
  IsLimitEnabled() {
    return this.m_enableLimit;
  }
  EnableLimit(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableLimit = flag;
  }
  GetLowerLimit() {
    return this.m_lowerTranslation;
  }
  GetUpperLimit() {
    return this.m_upperTranslation;
  }
  SetLimits(lower = 0, upper = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_lowerTranslation = lower;
    this.m_upperTranslation = upper;
  }
  IsMotorEnabled() {
    return this.m_enableMotor;
  }
  EnableMotor(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableMotor = flag;
  }
  SetMotorSpeed(speed = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_motorSpeed = speed;
  }
  GetMotorSpeed() {
    return this.m_motorSpeed;
  }
  SetMaxMotorForce(force = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_maxMotorForce = force;
  }
  GetMaxMotorForce() {
    return this.m_maxMotorForce;
  }
  GetMotorForce() {
    return this.m_motorImpulse;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tVec;
    this.m_localCenterA.SetV(bA.GetLocalCenter());
    this.m_localCenterB.SetV(bB.GetLocalCenter());
    const xf1 = bA.GetTransform();
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - this.m_localCenterA.x;
    let r1Y = this.m_localAnchor1.y - this.m_localCenterA.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - this.m_localCenterB.x;
    let r2Y = this.m_localAnchor2.y - this.m_localCenterB.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const dX = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
    const dY = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
    this.m_invMassA = bA.m_invMass;
    this.m_invMassB = bB.m_invMass;
    this.m_invIA = bA.m_invI;
    this.m_invIB = bB.m_invI;
    tMat = xf1.R;
    tVec = this.m_localXAxis1;
    this.m_axis.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_axis.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_a1 = (dX + r1X) * this.m_axis.y - (dY + r1Y) * this.m_axis.x;
    this.m_a2 = r2X * this.m_axis.y - r2Y * this.m_axis.x;
    this.m_motorMass = this.m_invMassA + this.m_invMassB + this.m_invIA * this.m_a1 * this.m_a1 + this.m_invIB * this.m_a2 * this.m_a2;
    this.m_motorMass = this.m_motorMass > b2Settings.b2_epsilon ? 1.0 / this.m_motorMass : 0.0;
    {
      tVec = this.m_localYAxis1;
      this.m_perp.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      this.m_perp.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      this.m_s1 = (dX + r1X) * this.m_perp.y - (dY + r1Y) * this.m_perp.x;
      this.m_s2 = r2X * this.m_perp.y - r2Y * this.m_perp.x;
      const m1 = this.m_invMassA;
      const m2 = this.m_invMassB;
      const i1 = this.m_invIA;
      const i2 = this.m_invIB;
      this.m_K.col1.x = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      this.m_K.col1.y = i1 * this.m_s1 * this.m_a1 + i2 * this.m_s2 * this.m_a2;
      this.m_K.col2.x = this.m_K.col1.y;
      this.m_K.col2.y = m1 + m2 + i1 * this.m_a1 * this.m_a1 + i2 * this.m_a2 * this.m_a2;
    }
    if (this.m_enableLimit) {
      const jointTransition = this.m_axis.x * dX + this.m_axis.y * dY;
      if (b2Math.Abs(this.m_upperTranslation - this.m_lowerTranslation) < 2.0 * b2Settings.b2_linearSlop) {
        this.m_limitState = b2Joint.e_equalLimits;
      } else if (jointTransition <= this.m_lowerTranslation) {
        if (this.m_limitState != b2Joint.e_atLowerLimit) {
          this.m_limitState = b2Joint.e_atLowerLimit;
          this.m_impulse.y = 0.0;
        }
      } else if (jointTransition >= this.m_upperTranslation) {
        if (this.m_limitState != b2Joint.e_atUpperLimit) {
          this.m_limitState = b2Joint.e_atUpperLimit;
          this.m_impulse.y = 0.0;
        }
      } else {
        this.m_limitState = b2Joint.e_inactiveLimit;
        this.m_impulse.y = 0.0;
      }
    } else {
      this.m_limitState = b2Joint.e_inactiveLimit;
    }
    if (this.m_enableMotor == false) {
      this.m_motorImpulse = 0.0;
    }
    if (step.warmStarting) {
      this.m_impulse.x *= step.dtRatio;
      this.m_impulse.y *= step.dtRatio;
      this.m_motorImpulse *= step.dtRatio;
      const PX = this.m_impulse.x * this.m_perp.x + (this.m_motorImpulse + this.m_impulse.y) * this.m_axis.x;
      const PY = this.m_impulse.x * this.m_perp.y + (this.m_motorImpulse + this.m_impulse.y) * this.m_axis.y;
      const L1 = this.m_impulse.x * this.m_s1 + (this.m_motorImpulse + this.m_impulse.y) * this.m_a1;
      const L2 = this.m_impulse.x * this.m_s2 + (this.m_motorImpulse + this.m_impulse.y) * this.m_a2;
      bA.m_linearVelocity.x -= this.m_invMassA * PX;
      bA.m_linearVelocity.y -= this.m_invMassA * PY;
      bA.m_angularVelocity -= this.m_invIA * L1;
      bB.m_linearVelocity.x += this.m_invMassB * PX;
      bB.m_linearVelocity.y += this.m_invMassB * PY;
      bB.m_angularVelocity += this.m_invIB * L2;
    } else {
      this.m_impulse.SetZero();
      this.m_motorImpulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const v1 = bA.m_linearVelocity;
    let w1 = bA.m_angularVelocity;
    const v2 = bB.m_linearVelocity;
    let w2 = bB.m_angularVelocity;
    let PX;
    let PY;
    let L1;
    let L2;
    if (this.m_enableMotor && this.m_limitState != b2Joint.e_equalLimits) {
      const Cdot = this.m_axis.x * (v2.x - v1.x) + this.m_axis.y * (v2.y - v1.y) + this.m_a2 * w2 - this.m_a1 * w1;
      let impulse = this.m_motorMass * (this.m_motorSpeed - Cdot);
      const oldImpulse = this.m_motorImpulse;
      const maxImpulse = step.dt * this.m_maxMotorForce;
      this.m_motorImpulse = b2Math.Clamp(this.m_motorImpulse + impulse, -maxImpulse, maxImpulse);
      impulse = this.m_motorImpulse - oldImpulse;
      PX = impulse * this.m_axis.x;
      PY = impulse * this.m_axis.y;
      L1 = impulse * this.m_a1;
      L2 = impulse * this.m_a2;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    }
    const Cdot1 = this.m_perp.x * (v2.x - v1.x) + this.m_perp.y * (v2.y - v1.y) + this.m_s2 * w2 - this.m_s1 * w1;
    if (this.m_enableLimit && this.m_limitState != b2Joint.e_inactiveLimit) {
      const Cdot2 = this.m_axis.x * (v2.x - v1.x) + this.m_axis.y * (v2.y - v1.y) + this.m_a2 * w2 - this.m_a1 * w1;
      const f1X = this.m_impulse.x;
      const f1Y = this.m_impulse.y;
      const df = this.m_K.Solve(this.m_tempImpulse, -Cdot1, -Cdot2);
      this.m_impulse.x += df.x;
      this.m_impulse.y += df.y;
      if (this.m_limitState == b2Joint.e_atLowerLimit) {
        this.m_impulse.y = b2Math.Max(this.m_impulse.y, 0.0);
      } else if (this.m_limitState == b2Joint.e_atUpperLimit) {
        this.m_impulse.y = b2Math.Min(this.m_impulse.y, 0.0);
      }
      const b = -Cdot1 - (this.m_impulse.y - f1Y) * this.m_K.col2.x;
      let f2r;
      if (this.m_K.col1.x != 0.0) {
        f2r = b / this.m_K.col1.x + f1X;
      } else {
        f2r = f1X;
      }
      this.m_impulse.x = f2r;
      df.x = this.m_impulse.x - f1X;
      df.y = this.m_impulse.y - f1Y;
      PX = df.x * this.m_perp.x + df.y * this.m_axis.x;
      PY = df.x * this.m_perp.y + df.y * this.m_axis.y;
      L1 = df.x * this.m_s1 + df.y * this.m_a1;
      L2 = df.x * this.m_s2 + df.y * this.m_a2;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    } else {
      let df2;
      if (this.m_K.col1.x != 0.0) {
        df2 = -Cdot1 / this.m_K.col1.x;
      } else {
        df2 = 0.0;
      }
      this.m_impulse.x += df2;
      PX = df2 * this.m_perp.x;
      PY = df2 * this.m_perp.y;
      L1 = df2 * this.m_s1;
      L2 = df2 * this.m_s2;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    }
    bA.m_angularVelocity = w1;
    bB.m_angularVelocity = w2;
  }
  SolvePositionConstraints(baumgarte) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const c1 = bA.m_sweep.c;
    let a1 = bA.m_sweep.a;
    const c2 = bB.m_sweep.c;
    let a2 = bB.m_sweep.a;
    let tVec;
    let m1;
    let m2;
    let i1;
    let i2;
    let linearError = 0.0;
    let active = false;
    let C2 = 0.0;
    const R1 = this.m_tempRotation1;
    const R2 = this.m_tempRotation2;
    R1.Set(a1);
    R2.Set(a2);
    let tMat = R1;
    let r1X = this.m_localAnchor1.x - this.m_localCenterA.x;
    let r1Y = this.m_localAnchor1.y - this.m_localCenterA.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = R2;
    let r2X = this.m_localAnchor2.x - this.m_localCenterB.x;
    let r2Y = this.m_localAnchor2.y - this.m_localCenterB.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const dX = c2.x + r2X - c1.x - r1X;
    const dY = c2.y + r2Y - c1.y - r1Y;
    if (this.m_enableLimit) {
      tMat = R1;
      tVec = this.m_localXAxis1;
      this.m_axis.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      this.m_axis.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      this.m_a1 = (dX + r1X) * this.m_axis.y - (dY + r1Y) * this.m_axis.x;
      this.m_a2 = r2X * this.m_axis.y - r2Y * this.m_axis.x;
      const translation = this.m_axis.x * dX + this.m_axis.y * dY;
      if (b2Math.Abs(this.m_upperTranslation - this.m_lowerTranslation) < 2.0 * b2Settings.b2_linearSlop) {
        C2 = b2Math.Clamp(translation, -b2Settings.b2_maxLinearCorrection, b2Settings.b2_maxLinearCorrection);
        linearError = b2Math.Abs(translation);
        active = true;
      } else if (translation <= this.m_lowerTranslation) {
        C2 = b2Math.Clamp(translation - this.m_lowerTranslation + b2Settings.b2_linearSlop, -b2Settings.b2_maxLinearCorrection, 0.0);
        linearError = this.m_lowerTranslation - translation;
        active = true;
      } else if (translation >= this.m_upperTranslation) {
        C2 = b2Math.Clamp(translation - this.m_upperTranslation - b2Settings.b2_linearSlop, 0.0, b2Settings.b2_maxLinearCorrection);
        linearError = translation - this.m_upperTranslation;
        active = true;
      }
    }
    tMat = R1;
    tVec = this.m_localYAxis1;
    this.m_perp.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_perp.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_s1 = (dX + r1X) * this.m_perp.y - (dY + r1Y) * this.m_perp.x;
    this.m_s2 = r2X * this.m_perp.y - r2Y * this.m_perp.x;
    const impulse = this.m_tempImpulse;
    const C1 = this.m_perp.x * dX + this.m_perp.y * dY;
    linearError = b2Math.Max(linearError, b2Math.Abs(C1));
    const angularError = 0.0;
    if (active) {
      m1 = this.m_invMassA;
      m2 = this.m_invMassB;
      i1 = this.m_invIA;
      i2 = this.m_invIB;
      this.m_K.col1.x = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      this.m_K.col1.y = i1 * this.m_s1 * this.m_a1 + i2 * this.m_s2 * this.m_a2;
      this.m_K.col2.x = this.m_K.col1.y;
      this.m_K.col2.y = m1 + m2 + i1 * this.m_a1 * this.m_a1 + i2 * this.m_a2 * this.m_a2;
      this.m_K.Solve(impulse, -C1, -C2);
    } else {
      m1 = this.m_invMassA;
      m2 = this.m_invMassB;
      i1 = this.m_invIA;
      i2 = this.m_invIB;
      const k11 = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      let impulse1;
      if (k11 != 0.0) {
        impulse1 = -C1 / k11;
      } else {
        impulse1 = 0.0;
      }
      impulse.x = impulse1;
      impulse.y = 0.0;
    }
    const PX = impulse.x * this.m_perp.x + impulse.y * this.m_axis.x;
    const PY = impulse.x * this.m_perp.y + impulse.y * this.m_axis.y;
    const L1 = impulse.x * this.m_s1 + impulse.y * this.m_a1;
    const L2 = impulse.x * this.m_s2 + impulse.y * this.m_a2;
    c1.x -= this.m_invMassA * PX;
    c1.y -= this.m_invMassA * PY;
    a1 -= this.m_invIA * L1;
    c2.x += this.m_invMassB * PX;
    c2.y += this.m_invMassB * PY;
    a2 += this.m_invIB * L2;
    bA.m_sweep.a = a1;
    bB.m_sweep.a = a2;
    bA.SynchronizeTransform();
    bB.SynchronizeTransform();
    return linearError <= b2Settings.b2_linearSlop && angularError <= b2Settings.b2_angularSlop;
  }
};

const b2LineJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.localAxisA = new b2Vec2();
    this.type = b2Joint.e_lineJoint;
    this.localAxisA.Set(1.0, 0.0);
    this.enableLimit = false;
    this.lowerTranslation = 0.0;
    this.upperTranslation = 0.0;
    this.enableMotor = false;
    this.maxMotorForce = 0.0;
    this.motorSpeed = 0.0;
  }
  Initialize(bA, bB, anchor, axis) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA = this.bodyA.GetLocalPoint(anchor);
    this.localAnchorB = this.bodyB.GetLocalPoint(anchor);
    this.localAxisA = this.bodyA.GetLocalVector(axis);
  }
};

const b2MouseJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.K = new b2Mat22();
    this.m_localAnchor = new b2Vec2();
    this.m_target = new b2Vec2();
    this.m_impulse = new b2Vec2();
    this.m_mass = new b2Mat22();
    this.m_C = new b2Vec2();
    this.m_target.SetV(def.target);
    const tX = this.m_target.x - this.m_bodyB.m_xf.position.x;
    const tY = this.m_target.y - this.m_bodyB.m_xf.position.y;
    const tMat = this.m_bodyB.m_xf.R;
    this.m_localAnchor.x = tX * tMat.col1.x + tY * tMat.col1.y;
    this.m_localAnchor.y = tX * tMat.col2.x + tY * tMat.col2.y;
    this.m_maxForce = def.maxForce;
    this.m_impulse.SetZero();
    this.m_frequencyHz = def.frequencyHz;
    this.m_dampingRatio = def.dampingRatio;
    this.m_beta = 0.0;
    this.m_gamma = 0.0;
  }
  GetAnchorA() {
    return this.m_target;
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse.x, inv_dt * this.m_impulse.y);
  }
  GetReactionTorque(inv_dt) {
    return 0.0;
  }
  GetTarget() {
    return this.m_target;
  }
  SetTarget(target) {
    if (this.m_bodyB.IsAwake() == false) {
      this.m_bodyB.SetAwake(true);
    }
    this.m_target.SetV(target);
  }
  GetMaxForce() {
    return this.m_maxForce;
  }
  SetMaxForce(maxForce = 0) {
    this.m_maxForce = maxForce;
  }
  GetFrequency() {
    return this.m_frequencyHz;
  }
  SetFrequency(hz = 0) {
    this.m_frequencyHz = hz;
  }
  GetDampingRatio() {
    return this.m_dampingRatio;
  }
  SetDampingRatio(ratio = 0) {
    this.m_dampingRatio = ratio;
  }
  InitVelocityConstraints(step) {
    const b = this.m_bodyB;
    const mass = b.m_mass;
    const omega = 2.0 * Math.PI * this.m_frequencyHz;
    const d = 2.0 * mass * this.m_dampingRatio * omega;
    const k = mass * omega * omega;
    this.m_gamma = step.dt * (d + step.dt * k);
    this.m_gamma = this.m_gamma != 0 ? 1 / this.m_gamma : 0.0;
    this.m_beta = step.dt * k * this.m_gamma;
    const tMat = b.m_xf.R;
    let rX = this.m_localAnchor.x - b.m_sweep.localCenter.x;
    let rY = this.m_localAnchor.y - b.m_sweep.localCenter.y;
    const tX = tMat.col1.x * rX + tMat.col2.x * rY;
    rY = tMat.col1.y * rX + tMat.col2.y * rY;
    rX = tX;
    const invMass = b.m_invMass;
    const invI = b.m_invI;
    this.K.col1.x = invMass + invI * rY * rY + this.m_gamma;
    this.K.col2.x = 0.0 + (-invI * rX * rY);
    this.K.col1.y = 0.0 + (-invI * rX * rY);
    this.K.col2.y = invMass + invI * rX * rX + this.m_gamma;
    this.K.GetInverse(this.m_mass);
    this.m_C.x = b.m_sweep.c.x + rX - this.m_target.x;
    this.m_C.y = b.m_sweep.c.y + rY - this.m_target.y;
    b.m_angularVelocity *= 0.98;
    this.m_impulse.x *= step.dtRatio;
    this.m_impulse.y *= step.dtRatio;
    b.m_linearVelocity.x += invMass * this.m_impulse.x;
    b.m_linearVelocity.y += invMass * this.m_impulse.y;
    b.m_angularVelocity += invI * (rX * this.m_impulse.y - rY * this.m_impulse.x);
  }
  SolveVelocityConstraints(step) {
    const b = this.m_bodyB;
    let tMat = b.m_xf.R;
    let rX = this.m_localAnchor.x - b.m_sweep.localCenter.x;
    let rY = this.m_localAnchor.y - b.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rX + tMat.col2.x * rY;
    rY = tMat.col1.y * rX + tMat.col2.y * rY;
    rX = tX;
    const CdotX = b.m_linearVelocity.x + (-b.m_angularVelocity * rY);
    const CdotY = b.m_linearVelocity.y + (b.m_angularVelocity * rX);
    tMat = this.m_mass;
    tX = CdotX + this.m_beta * this.m_C.x + this.m_gamma * this.m_impulse.x;
    const tY = CdotY + this.m_beta * this.m_C.y + this.m_gamma * this.m_impulse.y;
    let impulseX = -(tMat.col1.x * tX + tMat.col2.x * tY);
    let impulseY = -(tMat.col1.y * tX + tMat.col2.y * tY);
    const oldImpulseX = this.m_impulse.x;
    const oldImpulseY = this.m_impulse.y;
    this.m_impulse.x += impulseX;
    this.m_impulse.y += impulseY;
    const maxImpulse = step.dt * this.m_maxForce;
    if (this.m_impulse.LengthSquared() > maxImpulse * maxImpulse) {
      this.m_impulse.Multiply(maxImpulse / this.m_impulse.Length());
    }
    impulseX = this.m_impulse.x - oldImpulseX;
    impulseY = this.m_impulse.y - oldImpulseY;
    b.m_linearVelocity.x += b.m_invMass * impulseX;
    b.m_linearVelocity.y += b.m_invMass * impulseY;
    b.m_angularVelocity += b.m_invI * (rX * impulseY - rY * impulseX);
  }
  SolvePositionConstraints(baumgarte) {
    return true;
  }
};

const b2MouseJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.target = new b2Vec2();
    this.type = b2Joint.e_mouseJoint;
    this.maxForce = 0.0;
    this.frequencyHz = 5.0;
    this.dampingRatio = 0.7;
  }
};

const b2PrismaticJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_localXAxis1 = new b2Vec2();
    this.m_localYAxis1 = new b2Vec2();
    this.m_axis = new b2Vec2();
    this.m_perp = new b2Vec2();
    this.m_K = new b2Mat33();
    this.m_impulse = new b2Vec3();
    this.m_tempImpulse2 = new b2Vec2();
    this.m_tempImpulse3 = new b2Vec3();
    this.m_tempRotation1 = new b2Mat22();
    this.m_tempRotation2 = new b2Mat22();
    this.m_localAnchor1.SetV(def.localAnchorA);
    this.m_localAnchor2.SetV(def.localAnchorB);
    this.m_localXAxis1.SetV(def.localAxisA);
    this.m_localYAxis1.x = -this.m_localXAxis1.y;
    this.m_localYAxis1.y = this.m_localXAxis1.x;
    this.m_refAngle = def.referenceAngle;
    this.m_impulse.SetZero();
    this.m_motorMass = 0.0;
    this.m_motorImpulse = 0.0;
    this.m_lowerTranslation = def.lowerTranslation;
    this.m_upperTranslation = def.upperTranslation;
    this.m_maxMotorForce = def.maxMotorForce;
    this.m_motorSpeed = def.motorSpeed;
    this.m_enableLimit = def.enableLimit;
    this.m_enableMotor = def.enableMotor;
    this.m_limitState = b2Joint.e_inactiveLimit;
    this.m_axis.SetZero();
    this.m_perp.SetZero();
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * (this.m_impulse.x * this.m_perp.x + (this.m_motorImpulse + this.m_impulse.z) * this.m_axis.x), inv_dt * (this.m_impulse.x * this.m_perp.y + (this.m_motorImpulse + this.m_impulse.z) * this.m_axis.y));
  }
  GetReactionTorque(inv_dt = 0) {
    return inv_dt * this.m_impulse.y;
  }
  GetJointTranslation() {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let tVec = this.m_localAnchor1;
    const p1X = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + bA.m_xf.position.x;
    const p1Y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + bA.m_xf.position.y;
    tVec = this.m_localXAxis1;
    const axisX = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    const axisY = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    tMat = bB.m_xf.R;
    tVec = this.m_localAnchor2;
    const p2X = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y + bB.m_xf.position.x;
    const p2Y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y + bB.m_xf.position.y;
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    const translation = axisX * dX + axisY * dY;
    return translation;
  }
  GetJointSpeed() {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const p1X = bA.m_sweep.c.x + r1X;
    const p1Y = bA.m_sweep.c.y + r1Y;
    const p2X = bB.m_sweep.c.x + r2X;
    const p2Y = bB.m_sweep.c.y + r2Y;
    const dX = p2X - p1X;
    const dY = p2Y - p1Y;
    const axis = bA.GetWorldVector(this.m_localXAxis1);
    const v1 = bA.m_linearVelocity;
    const v2 = bB.m_linearVelocity;
    const w1 = bA.m_angularVelocity;
    const w2 = bB.m_angularVelocity;
    const speed = (dX * (-w1 * axis.y) + dY * (w1 * axis.x)) + (axis.x * (((v2.x + (-w2 * r2Y)) - v1.x) - (-w1 * r1Y)) + axis.y * (((v2.y + (w2 * r2X)) - v1.y) - (w1 * r1X)));
    return speed;
  }
  IsLimitEnabled() {
    return this.m_enableLimit;
  }
  EnableLimit(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableLimit = flag;
  }
  GetLowerLimit() {
    return this.m_lowerTranslation;
  }
  GetUpperLimit() {
    return this.m_upperTranslation;
  }
  SetLimits(lower = 0, upper = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_lowerTranslation = lower;
    this.m_upperTranslation = upper;
  }
  IsMotorEnabled() {
    return this.m_enableMotor;
  }
  EnableMotor(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableMotor = flag;
  }
  SetMotorSpeed(speed = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_motorSpeed = speed;
  }
  GetMotorSpeed() {
    return this.m_motorSpeed;
  }
  SetMaxMotorForce(force = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_maxMotorForce = force;
  }
  GetMotorForce() {
    return this.m_motorImpulse;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tVec;
    this.m_localCenterA.SetV(bA.GetLocalCenter());
    this.m_localCenterB.SetV(bB.GetLocalCenter());
    const xf1 = bA.GetTransform();
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - this.m_localCenterA.x;
    let r1Y = this.m_localAnchor1.y - this.m_localCenterA.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - this.m_localCenterB.x;
    let r2Y = this.m_localAnchor2.y - this.m_localCenterB.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const dX = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
    const dY = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
    this.m_invMassA = bA.m_invMass;
    this.m_invMassB = bB.m_invMass;
    this.m_invIA = bA.m_invI;
    this.m_invIB = bB.m_invI;
    tMat = xf1.R;
    tVec = this.m_localXAxis1;
    this.m_axis.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_axis.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_a1 = (dX + r1X) * this.m_axis.y - (dY + r1Y) * this.m_axis.x;
    this.m_a2 = r2X * this.m_axis.y - r2Y * this.m_axis.x;
    this.m_motorMass = this.m_invMassA + this.m_invMassB + this.m_invIA * this.m_a1 * this.m_a1 + this.m_invIB * this.m_a2 * this.m_a2;
    if (this.m_motorMass > b2Settings.b2_epsilon) this.m_motorMass = 1.0 / this.m_motorMass;
    {
      tVec = this.m_localYAxis1;
      this.m_perp.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      this.m_perp.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      this.m_s1 = (dX + r1X) * this.m_perp.y - (dY + r1Y) * this.m_perp.x;
      this.m_s2 = r2X * this.m_perp.y - r2Y * this.m_perp.x;
      const m1 = this.m_invMassA;
      const m2 = this.m_invMassB;
      const i1 = this.m_invIA;
      const i2 = this.m_invIB;
      this.m_K.col1.x = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      this.m_K.col1.y = i1 * this.m_s1 + i2 * this.m_s2;
      this.m_K.col1.z = i1 * this.m_s1 * this.m_a1 + i2 * this.m_s2 * this.m_a2;
      this.m_K.col2.x = this.m_K.col1.y;
      this.m_K.col2.y = i1 + i2;
      this.m_K.col2.z = i1 * this.m_a1 + i2 * this.m_a2;
      this.m_K.col3.x = this.m_K.col1.z;
      this.m_K.col3.y = this.m_K.col2.z;
      this.m_K.col3.z = m1 + m2 + i1 * this.m_a1 * this.m_a1 + i2 * this.m_a2 * this.m_a2;
    }
    if (this.m_enableLimit) {
      const jointTransition = this.m_axis.x * dX + this.m_axis.y * dY;
      if (b2Math.Abs(this.m_upperTranslation - this.m_lowerTranslation) < 2.0 * b2Settings.b2_linearSlop) {
        this.m_limitState = b2Joint.e_equalLimits;
      } else if (jointTransition <= this.m_lowerTranslation) {
        if (this.m_limitState != b2Joint.e_atLowerLimit) {
          this.m_limitState = b2Joint.e_atLowerLimit;
          this.m_impulse.z = 0.0;
        }
      } else if (jointTransition >= this.m_upperTranslation) {
        if (this.m_limitState != b2Joint.e_atUpperLimit) {
          this.m_limitState = b2Joint.e_atUpperLimit;
          this.m_impulse.z = 0.0;
        }
      } else {
        this.m_limitState = b2Joint.e_inactiveLimit;
        this.m_impulse.z = 0.0;
      }
    } else {
      this.m_limitState = b2Joint.e_inactiveLimit;
      this.m_impulse.z = 0.0;
    }
    if (this.m_enableMotor == false) {
      this.m_motorImpulse = 0.0;
    }
    if (step.warmStarting) {
      this.m_impulse.x *= step.dtRatio;
      this.m_impulse.y *= step.dtRatio;
      this.m_impulse.z *= step.dtRatio;
      this.m_motorImpulse *= step.dtRatio;
      const PX = this.m_impulse.x * this.m_perp.x + (this.m_motorImpulse + this.m_impulse.z) * this.m_axis.x;
      const PY = this.m_impulse.x * this.m_perp.y + (this.m_motorImpulse + this.m_impulse.z) * this.m_axis.y;
      const L1 = this.m_impulse.x * this.m_s1 + this.m_impulse.y + (this.m_motorImpulse + this.m_impulse.z) * this.m_a1;
      const L2 = this.m_impulse.x * this.m_s2 + this.m_impulse.y + (this.m_motorImpulse + this.m_impulse.z) * this.m_a2;
      bA.m_linearVelocity.x -= this.m_invMassA * PX;
      bA.m_linearVelocity.y -= this.m_invMassA * PY;
      bA.m_angularVelocity -= this.m_invIA * L1;
      bB.m_linearVelocity.x += this.m_invMassB * PX;
      bB.m_linearVelocity.y += this.m_invMassB * PY;
      bB.m_angularVelocity += this.m_invIB * L2;
    } else {
      this.m_impulse.SetZero();
      this.m_motorImpulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const v1 = bA.m_linearVelocity;
    let w1 = bA.m_angularVelocity;
    const v2 = bB.m_linearVelocity;
    let w2 = bB.m_angularVelocity;
    let PX;
    let PY;
    let L1;
    let L2;
    if (this.m_enableMotor && this.m_limitState != b2Joint.e_equalLimits) {
      const Cdot = this.m_axis.x * (v2.x - v1.x) + this.m_axis.y * (v2.y - v1.y) + this.m_a2 * w2 - this.m_a1 * w1;
      let impulse = this.m_motorMass * (this.m_motorSpeed - Cdot);
      const oldImpulse = this.m_motorImpulse;
      const maxImpulse = step.dt * this.m_maxMotorForce;
      this.m_motorImpulse = b2Math.Clamp(this.m_motorImpulse + impulse, -maxImpulse, maxImpulse);
      impulse = this.m_motorImpulse - oldImpulse;
      PX = impulse * this.m_axis.x;
      PY = impulse * this.m_axis.y;
      L1 = impulse * this.m_a1;
      L2 = impulse * this.m_a2;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    }
    const Cdot1X = this.m_perp.x * (v2.x - v1.x) + this.m_perp.y * (v2.y - v1.y) + this.m_s2 * w2 - this.m_s1 * w1;
    const Cdot1Y = w2 - w1;
    if (this.m_enableLimit && this.m_limitState != b2Joint.e_inactiveLimit) {
      const Cdot2 = this.m_axis.x * (v2.x - v1.x) + this.m_axis.y * (v2.y - v1.y) + this.m_a2 * w2 - this.m_a1 * w1;
      const f1X = this.m_impulse.x;
      const f1Y = this.m_impulse.y;
      const f1Z = this.m_impulse.z;
      const df = this.m_K.Solve33(this.m_tempImpulse3, -Cdot1X, -Cdot1Y, -Cdot2);
      this.m_impulse.x += df.x;
      this.m_impulse.y += df.y;
      this.m_impulse.z += df.z;
      if (this.m_limitState == b2Joint.e_atLowerLimit) {
        this.m_impulse.z = b2Math.Max(this.m_impulse.z, 0.0);
      } else if (this.m_limitState == b2Joint.e_atUpperLimit) {
        this.m_impulse.z = b2Math.Min(this.m_impulse.z, 0.0);
      }
      const bX = -Cdot1X - (this.m_impulse.z - f1Z) * this.m_K.col3.x;
      const bY = -Cdot1Y - (this.m_impulse.z - f1Z) * this.m_K.col3.y;
      const f2r = this.m_K.Solve22(this.m_tempImpulse2, bX, bY);
      f2r.x += f1X;
      f2r.y += f1Y;
      this.m_impulse.x = f2r.x;
      this.m_impulse.y = f2r.y;
      df.x = this.m_impulse.x - f1X;
      df.y = this.m_impulse.y - f1Y;
      df.z = this.m_impulse.z - f1Z;
      PX = df.x * this.m_perp.x + df.z * this.m_axis.x;
      PY = df.x * this.m_perp.y + df.z * this.m_axis.y;
      L1 = df.x * this.m_s1 + df.y + df.z * this.m_a1;
      L2 = df.x * this.m_s2 + df.y + df.z * this.m_a2;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    } else {
      const df2 = this.m_K.Solve22(this.m_tempImpulse2, -Cdot1X, -Cdot1Y);
      this.m_impulse.x += df2.x;
      this.m_impulse.y += df2.y;
      PX = df2.x * this.m_perp.x;
      PY = df2.x * this.m_perp.y;
      L1 = df2.x * this.m_s1 + df2.y;
      L2 = df2.x * this.m_s2 + df2.y;
      v1.x -= this.m_invMassA * PX;
      v1.y -= this.m_invMassA * PY;
      w1 -= this.m_invIA * L1;
      v2.x += this.m_invMassB * PX;
      v2.y += this.m_invMassB * PY;
      w2 += this.m_invIB * L2;
    }
    bA.m_angularVelocity = w1;
    bB.m_angularVelocity = w2;
  }
  SolvePositionConstraints(baumgarte) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const c1 = bA.m_sweep.c;
    let a1 = bA.m_sweep.a;
    const c2 = bB.m_sweep.c;
    let a2 = bB.m_sweep.a;
    let tVec;
    let m1;
    let m2;
    let i1;
    let i2;
    let linearError = 0.0;
    let active = false;
    let C2 = 0.0;
    const R1 = this.m_tempRotation1;
    const R2 = this.m_tempRotation2;
    R1.Set(a1);
    R2.Set(a2);
    let tMat = R1;
    let r1X = this.m_localAnchor1.x - this.m_localCenterA.x;
    let r1Y = this.m_localAnchor1.y - this.m_localCenterA.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = R2;
    let r2X = this.m_localAnchor2.x - this.m_localCenterB.x;
    let r2Y = this.m_localAnchor2.y - this.m_localCenterB.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const dX = c2.x + r2X - c1.x - r1X;
    const dY = c2.y + r2Y - c1.y - r1Y;
    if (this.m_enableLimit) {
      tMat = R1;
      tVec = this.m_localXAxis1;
      this.m_axis.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
      this.m_axis.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
      this.m_a1 = (dX + r1X) * this.m_axis.y - (dY + r1Y) * this.m_axis.x;
      this.m_a2 = r2X * this.m_axis.y - r2Y * this.m_axis.x;
      const translation = this.m_axis.x * dX + this.m_axis.y * dY;
      if (b2Math.Abs(this.m_upperTranslation - this.m_lowerTranslation) < 2.0 * b2Settings.b2_linearSlop) {
        C2 = b2Math.Clamp(translation, -b2Settings.b2_maxLinearCorrection, b2Settings.b2_maxLinearCorrection);
        linearError = b2Math.Abs(translation);
        active = true;
      } else if (translation <= this.m_lowerTranslation) {
        C2 = b2Math.Clamp(translation - this.m_lowerTranslation + b2Settings.b2_linearSlop, -b2Settings.b2_maxLinearCorrection, 0.0);
        linearError = this.m_lowerTranslation - translation;
        active = true;
      } else if (translation >= this.m_upperTranslation) {
        C2 = b2Math.Clamp(translation - this.m_upperTranslation - b2Settings.b2_linearSlop, 0.0, b2Settings.b2_maxLinearCorrection);
        linearError = translation - this.m_upperTranslation;
        active = true;
      }
    }
    tMat = R1;
    tVec = this.m_localYAxis1;
    this.m_perp.x = tMat.col1.x * tVec.x + tMat.col2.x * tVec.y;
    this.m_perp.y = tMat.col1.y * tVec.x + tMat.col2.y * tVec.y;
    this.m_s1 = (dX + r1X) * this.m_perp.y - (dY + r1Y) * this.m_perp.x;
    this.m_s2 = r2X * this.m_perp.y - r2Y * this.m_perp.x;
    const impulse = this.m_tempImpulse3;
    const C1X = this.m_perp.x * dX + this.m_perp.y * dY;
    const C1Y = a2 - a1 - this.m_refAngle;
    linearError = b2Math.Max(linearError, b2Math.Abs(C1X));
    const angularError = b2Math.Abs(C1Y);
    if (active) {
      m1 = this.m_invMassA;
      m2 = this.m_invMassB;
      i1 = this.m_invIA;
      i2 = this.m_invIB;
      this.m_K.col1.x = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      this.m_K.col1.y = i1 * this.m_s1 + i2 * this.m_s2;
      this.m_K.col1.z = i1 * this.m_s1 * this.m_a1 + i2 * this.m_s2 * this.m_a2;
      this.m_K.col2.x = this.m_K.col1.y;
      this.m_K.col2.y = i1 + i2;
      this.m_K.col2.z = i1 * this.m_a1 + i2 * this.m_a2;
      this.m_K.col3.x = this.m_K.col1.z;
      this.m_K.col3.y = this.m_K.col2.z;
      this.m_K.col3.z = m1 + m2 + i1 * this.m_a1 * this.m_a1 + i2 * this.m_a2 * this.m_a2;
      this.m_K.Solve33(impulse, -C1X, -C1Y, -C2);
    } else {
      m1 = this.m_invMassA;
      m2 = this.m_invMassB;
      i1 = this.m_invIA;
      i2 = this.m_invIB;
      const k11 = m1 + m2 + i1 * this.m_s1 * this.m_s1 + i2 * this.m_s2 * this.m_s2;
      const k12 = i1 * this.m_s1 + i2 * this.m_s2;
      const k22 = i1 + i2;
      this.m_K.col1.Set(k11, k12, 0.0);
      this.m_K.col2.Set(k12, k22, 0.0);
      const impulse1 = this.m_K.Solve22(this.m_tempImpulse2, -C1X, -C1Y);
      impulse.x = impulse1.x;
      impulse.y = impulse1.y;
      impulse.z = 0.0;
    }
    const PX = impulse.x * this.m_perp.x + impulse.z * this.m_axis.x;
    const PY = impulse.x * this.m_perp.y + impulse.z * this.m_axis.y;
    const L1 = impulse.x * this.m_s1 + impulse.y + impulse.z * this.m_a1;
    const L2 = impulse.x * this.m_s2 + impulse.y + impulse.z * this.m_a2;
    c1.x -= this.m_invMassA * PX;
    c1.y -= this.m_invMassA * PY;
    a1 -= this.m_invIA * L1;
    c2.x += this.m_invMassB * PX;
    c2.y += this.m_invMassB * PY;
    a2 += this.m_invIB * L2;
    bA.m_sweep.a = a1;
    bB.m_sweep.a = a2;
    bA.SynchronizeTransform();
    bB.SynchronizeTransform();
    return linearError <= b2Settings.b2_linearSlop && angularError <= b2Settings.b2_angularSlop;
  }
};

const b2PrismaticJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.localAxisA = new b2Vec2();
    this.type = b2Joint.e_prismaticJoint;
    this.localAxisA.Set(1.0, 0.0);
    this.referenceAngle = 0.0;
    this.enableLimit = false;
    this.lowerTranslation = 0.0;
    this.upperTranslation = 0.0;
    this.enableMotor = false;
    this.maxMotorForce = 0.0;
    this.motorSpeed = 0.0;
  }
  Initialize(bA, bB, anchor, axis) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA = this.bodyA.GetLocalPoint(anchor);
    this.localAnchorB = this.bodyB.GetLocalPoint(anchor);
    this.localAxisA = this.bodyA.GetLocalVector(axis);
    this.referenceAngle = this.bodyB.GetAngle() - this.bodyA.GetAngle();
  }
};

const b2PulleyJoint = class extends b2Joint {
  static b2_minPulleyLength = 2.0;
  constructor(def) {
    super(def);
    this.m_groundAnchor1 = new b2Vec2();
    this.m_groundAnchor2 = new b2Vec2();
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_u1 = new b2Vec2();
    this.m_u2 = new b2Vec2();
    this.m_ground = this.m_bodyA.m_world.m_groundBody;
    this.m_groundAnchor1.x = def.groundAnchorA.x - this.m_ground.m_xf.position.x;
    this.m_groundAnchor1.y = def.groundAnchorA.y - this.m_ground.m_xf.position.y;
    this.m_groundAnchor2.x = def.groundAnchorB.x - this.m_ground.m_xf.position.x;
    this.m_groundAnchor2.y = def.groundAnchorB.y - this.m_ground.m_xf.position.y;
    this.m_localAnchor1.SetV(def.localAnchorA);
    this.m_localAnchor2.SetV(def.localAnchorB);
    this.m_ratio = def.ratio;
    this.m_constant = def.lengthA + this.m_ratio * def.lengthB;
    this.m_maxLength1 = b2Math.Min(def.maxLengthA, this.m_constant - this.m_ratio * b2PulleyJoint.b2_minPulleyLength);
    this.m_maxLength2 = b2Math.Min(def.maxLengthB, (this.m_constant - b2PulleyJoint.b2_minPulleyLength) / this.m_ratio);
    this.m_impulse = 0.0;
    this.m_limitImpulse1 = 0.0;
    this.m_limitImpulse2 = 0.0;
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse * this.m_u2.x, inv_dt * this.m_impulse * this.m_u2.y);
  }
  GetReactionTorque(inv_dt) {
    return 0.0;
  }
  GetGroundAnchorA() {
    const a = this.m_ground.m_xf.position.Copy();
    a.Add(this.m_groundAnchor1);
    return a;
  }
  GetGroundAnchorB() {
    const a = this.m_ground.m_xf.position.Copy();
    a.Add(this.m_groundAnchor2);
    return a;
  }
  GetLength1() {
    const p = this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
    const sX = this.m_ground.m_xf.position.x + this.m_groundAnchor1.x;
    const sY = this.m_ground.m_xf.position.y + this.m_groundAnchor1.y;
    const dX = p.x - sX;
    const dY = p.y - sY;
    return Math.sqrt(dX * dX + dY * dY);
  }
  GetLength2() {
    const p = this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
    const sX = this.m_ground.m_xf.position.x + this.m_groundAnchor2.x;
    const sY = this.m_ground.m_xf.position.y + this.m_groundAnchor2.y;
    const dX = p.x - sX;
    const dY = p.y - sY;
    return Math.sqrt(dX * dX + dY * dY);
  }
  GetRatio() {
    return this.m_ratio;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const p1X = bA.m_sweep.c.x + r1X;
    const p1Y = bA.m_sweep.c.y + r1Y;
    const p2X = bB.m_sweep.c.x + r2X;
    const p2Y = bB.m_sweep.c.y + r2Y;
    const s1X = this.m_ground.m_xf.position.x + this.m_groundAnchor1.x;
    const s1Y = this.m_ground.m_xf.position.y + this.m_groundAnchor1.y;
    const s2X = this.m_ground.m_xf.position.x + this.m_groundAnchor2.x;
    const s2Y = this.m_ground.m_xf.position.y + this.m_groundAnchor2.y;
    this.m_u1.Set(p1X - s1X, p1Y - s1Y);
    this.m_u2.Set(p2X - s2X, p2Y - s2Y);
    const length1 = this.m_u1.Length();
    const length2 = this.m_u2.Length();
    if (length1 > b2Settings.b2_linearSlop) {
      this.m_u1.Multiply(1.0 / length1);
    } else {
      this.m_u1.SetZero();
    }
    if (length2 > b2Settings.b2_linearSlop) {
      this.m_u2.Multiply(1.0 / length2);
    } else {
      this.m_u2.SetZero();
    }
    const C = this.m_constant - length1 - this.m_ratio * length2;
    if (C > 0.0) {
      this.m_state = b2Joint.e_inactiveLimit;
      this.m_impulse = 0.0;
    } else {
      this.m_state = b2Joint.e_atUpperLimit;
    }
    if (length1 < this.m_maxLength1) {
      this.m_limitState1 = b2Joint.e_inactiveLimit;
      this.m_limitImpulse1 = 0.0;
    } else {
      this.m_limitState1 = b2Joint.e_atUpperLimit;
    }
    if (length2 < this.m_maxLength2) {
      this.m_limitState2 = b2Joint.e_inactiveLimit;
      this.m_limitImpulse2 = 0.0;
    } else {
      this.m_limitState2 = b2Joint.e_atUpperLimit;
    }
    const cr1u1 = r1X * this.m_u1.y - r1Y * this.m_u1.x;
    const cr2u2 = r2X * this.m_u2.y - r2Y * this.m_u2.x;
    this.m_limitMass1 = bA.m_invMass + bA.m_invI * cr1u1 * cr1u1;
    this.m_limitMass2 = bB.m_invMass + bB.m_invI * cr2u2 * cr2u2;
    this.m_pulleyMass = this.m_limitMass1 + this.m_ratio * this.m_ratio * this.m_limitMass2;
    this.m_limitMass1 = 1.0 / this.m_limitMass1;
    this.m_limitMass2 = 1.0 / this.m_limitMass2;
    this.m_pulleyMass = 1.0 / this.m_pulleyMass;
    if (step.warmStarting) {
      this.m_impulse *= step.dtRatio;
      this.m_limitImpulse1 *= step.dtRatio;
      this.m_limitImpulse2 *= step.dtRatio;
      const P1X = (-this.m_impulse - this.m_limitImpulse1) * this.m_u1.x;
      const P1Y = (-this.m_impulse - this.m_limitImpulse1) * this.m_u1.y;
      const P2X = ((-this.m_ratio * this.m_impulse) - this.m_limitImpulse2) * this.m_u2.x;
      const P2Y = ((-this.m_ratio * this.m_impulse) - this.m_limitImpulse2) * this.m_u2.y;
      bA.m_linearVelocity.x += bA.m_invMass * P1X;
      bA.m_linearVelocity.y += bA.m_invMass * P1Y;
      bA.m_angularVelocity += bA.m_invI * (r1X * P1Y - r1Y * P1X);
      bB.m_linearVelocity.x += bB.m_invMass * P2X;
      bB.m_linearVelocity.y += bB.m_invMass * P2Y;
      bB.m_angularVelocity += bB.m_invI * (r2X * P2Y - r2Y * P2X);
    } else {
      this.m_impulse = 0.0;
      this.m_limitImpulse1 = 0.0;
      this.m_limitImpulse2 = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    let v1X;
    let v1Y;
    let v2X;
    let v2Y;
    let P1X;
    let P1Y;
    let P2X;
    let P2Y;
    let Cdot;
    let impulse;
    let oldImpulse;
    if (this.m_state == b2Joint.e_atUpperLimit) {
      v1X = bA.m_linearVelocity.x + (-bA.m_angularVelocity * r1Y);
      v1Y = bA.m_linearVelocity.y + (bA.m_angularVelocity * r1X);
      v2X = bB.m_linearVelocity.x + (-bB.m_angularVelocity * r2Y);
      v2Y = bB.m_linearVelocity.y + (bB.m_angularVelocity * r2X);
      Cdot = -(this.m_u1.x * v1X + this.m_u1.y * v1Y) - this.m_ratio * (this.m_u2.x * v2X + this.m_u2.y * v2Y);
      impulse = this.m_pulleyMass * -Cdot;
      oldImpulse = this.m_impulse;
      this.m_impulse = b2Math.Max(0.0, this.m_impulse + impulse);
      impulse = this.m_impulse - oldImpulse;
      P1X = -impulse * this.m_u1.x;
      P1Y = -impulse * this.m_u1.y;
      P2X = -this.m_ratio * impulse * this.m_u2.x;
      P2Y = -this.m_ratio * impulse * this.m_u2.y;
      bA.m_linearVelocity.x += bA.m_invMass * P1X;
      bA.m_linearVelocity.y += bA.m_invMass * P1Y;
      bA.m_angularVelocity += bA.m_invI * (r1X * P1Y - r1Y * P1X);
      bB.m_linearVelocity.x += bB.m_invMass * P2X;
      bB.m_linearVelocity.y += bB.m_invMass * P2Y;
      bB.m_angularVelocity += bB.m_invI * (r2X * P2Y - r2Y * P2X);
    }
    if (this.m_limitState1 == b2Joint.e_atUpperLimit) {
      v1X = bA.m_linearVelocity.x + (-bA.m_angularVelocity * r1Y);
      v1Y = bA.m_linearVelocity.y + (bA.m_angularVelocity * r1X);
      Cdot = -(this.m_u1.x * v1X + this.m_u1.y * v1Y);
      impulse = -this.m_limitMass1 * Cdot;
      oldImpulse = this.m_limitImpulse1;
      this.m_limitImpulse1 = b2Math.Max(0.0, this.m_limitImpulse1 + impulse);
      impulse = this.m_limitImpulse1 - oldImpulse;
      P1X = -impulse * this.m_u1.x;
      P1Y = -impulse * this.m_u1.y;
      bA.m_linearVelocity.x += bA.m_invMass * P1X;
      bA.m_linearVelocity.y += bA.m_invMass * P1Y;
      bA.m_angularVelocity += bA.m_invI * (r1X * P1Y - r1Y * P1X);
    }
    if (this.m_limitState2 == b2Joint.e_atUpperLimit) {
      v2X = bB.m_linearVelocity.x + (-bB.m_angularVelocity * r2Y);
      v2Y = bB.m_linearVelocity.y + (bB.m_angularVelocity * r2X);
      Cdot = -(this.m_u2.x * v2X + this.m_u2.y * v2Y);
      impulse = -this.m_limitMass2 * Cdot;
      oldImpulse = this.m_limitImpulse2;
      this.m_limitImpulse2 = b2Math.Max(0.0, this.m_limitImpulse2 + impulse);
      impulse = this.m_limitImpulse2 - oldImpulse;
      P2X = -impulse * this.m_u2.x;
      P2Y = -impulse * this.m_u2.y;
      bB.m_linearVelocity.x += bB.m_invMass * P2X;
      bB.m_linearVelocity.y += bB.m_invMass * P2Y;
      bB.m_angularVelocity += bB.m_invI * (r2X * P2Y - r2Y * P2X);
    }
  }
  SolvePositionConstraints(baumgarte) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat;
    const s1X = this.m_ground.m_xf.position.x + this.m_groundAnchor1.x;
    const s1Y = this.m_ground.m_xf.position.y + this.m_groundAnchor1.y;
    const s2X = this.m_ground.m_xf.position.x + this.m_groundAnchor2.x;
    const s2Y = this.m_ground.m_xf.position.y + this.m_groundAnchor2.y;
    let r1X;
    let r1Y;
    let r2X;
    let r2Y;
    let p1X;
    let p1Y;
    let p2X;
    let p2Y;
    let length1;
    let length2;
    let C;
    let impulse;
    let tX;
    let linearError = 0.0;
    if (this.m_state == b2Joint.e_atUpperLimit) {
      tMat = bA.m_xf.R;
      r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
      r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
      r1X = tX;
      tMat = bB.m_xf.R;
      r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
      r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
      r2X = tX;
      p1X = bA.m_sweep.c.x + r1X;
      p1Y = bA.m_sweep.c.y + r1Y;
      p2X = bB.m_sweep.c.x + r2X;
      p2Y = bB.m_sweep.c.y + r2Y;
      this.m_u1.Set(p1X - s1X, p1Y - s1Y);
      this.m_u2.Set(p2X - s2X, p2Y - s2Y);
      length1 = this.m_u1.Length();
      length2 = this.m_u2.Length();
      if (length1 > b2Settings.b2_linearSlop) {
        this.m_u1.Multiply(1.0 / length1);
      } else {
        this.m_u1.SetZero();
      }
      if (length2 > b2Settings.b2_linearSlop) {
        this.m_u2.Multiply(1.0 / length2);
      } else {
        this.m_u2.SetZero();
      }
      C = this.m_constant - length1 - this.m_ratio * length2;
      linearError = b2Math.Max(linearError, -C);
      C = b2Math.Clamp(C + b2Settings.b2_linearSlop, -b2Settings.b2_maxLinearCorrection, 0.0);
      impulse = -this.m_pulleyMass * C;
      p1X = -impulse * this.m_u1.x;
      p1Y = -impulse * this.m_u1.y;
      p2X = -this.m_ratio * impulse * this.m_u2.x;
      p2Y = -this.m_ratio * impulse * this.m_u2.y;
      bA.m_sweep.c.x += bA.m_invMass * p1X;
      bA.m_sweep.c.y += bA.m_invMass * p1Y;
      bA.m_sweep.a += bA.m_invI * (r1X * p1Y - r1Y * p1X);
      bB.m_sweep.c.x += bB.m_invMass * p2X;
      bB.m_sweep.c.y += bB.m_invMass * p2Y;
      bB.m_sweep.a += bB.m_invI * (r2X * p2Y - r2Y * p2X);
      bA.SynchronizeTransform();
      bB.SynchronizeTransform();
    }
    if (this.m_limitState1 == b2Joint.e_atUpperLimit) {
      tMat = bA.m_xf.R;
      r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
      r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
      r1X = tX;
      p1X = bA.m_sweep.c.x + r1X;
      p1Y = bA.m_sweep.c.y + r1Y;
      this.m_u1.Set(p1X - s1X, p1Y - s1Y);
      length1 = this.m_u1.Length();
      if (length1 > b2Settings.b2_linearSlop) {
        this.m_u1.x *= 1.0 / length1;
        this.m_u1.y *= 1.0 / length1;
      } else {
        this.m_u1.SetZero();
      }
      C = this.m_maxLength1 - length1;
      linearError = b2Math.Max(linearError, -C);
      C = b2Math.Clamp(C + b2Settings.b2_linearSlop, -b2Settings.b2_maxLinearCorrection, 0.0);
      impulse = -this.m_limitMass1 * C;
      p1X = -impulse * this.m_u1.x;
      p1Y = -impulse * this.m_u1.y;
      bA.m_sweep.c.x += bA.m_invMass * p1X;
      bA.m_sweep.c.y += bA.m_invMass * p1Y;
      bA.m_sweep.a += bA.m_invI * (r1X * p1Y - r1Y * p1X);
      bA.SynchronizeTransform();
    }
    if (this.m_limitState2 == b2Joint.e_atUpperLimit) {
      tMat = bB.m_xf.R;
      r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
      r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
      r2X = tX;
      p2X = bB.m_sweep.c.x + r2X;
      p2Y = bB.m_sweep.c.y + r2Y;
      this.m_u2.Set(p2X - s2X, p2Y - s2Y);
      length2 = this.m_u2.Length();
      if (length2 > b2Settings.b2_linearSlop) {
        this.m_u2.x *= 1.0 / length2;
        this.m_u2.y *= 1.0 / length2;
      } else {
        this.m_u2.SetZero();
      }
      C = this.m_maxLength2 - length2;
      linearError = b2Math.Max(linearError, -C);
      C = b2Math.Clamp(C + b2Settings.b2_linearSlop, -b2Settings.b2_maxLinearCorrection, 0.0);
      impulse = -this.m_limitMass2 * C;
      p2X = -impulse * this.m_u2.x;
      p2Y = -impulse * this.m_u2.y;
      bB.m_sweep.c.x += bB.m_invMass * p2X;
      bB.m_sweep.c.y += bB.m_invMass * p2Y;
      bB.m_sweep.a += bB.m_invI * (r2X * p2Y - r2Y * p2X);
      bB.SynchronizeTransform();
    }
    return linearError < b2Settings.b2_linearSlop;
  }
};

const b2PulleyJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.groundAnchorA = new b2Vec2();
    this.groundAnchorB = new b2Vec2();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.type = b2Joint.e_pulleyJoint;
    this.groundAnchorA.Set(-1.0, 1.0);
    this.groundAnchorB.Set(1.0, 1.0);
    this.localAnchorA.Set(-1.0, 0.0);
    this.localAnchorB.Set(1.0, 0.0);
    this.lengthA = 0.0;
    this.maxLengthA = 0.0;
    this.lengthB = 0.0;
    this.maxLengthB = 0.0;
    this.ratio = 1.0;
    this.collideConnected = true;
  }
  Initialize(bA, bB, gaA, gaB, anchorA, anchorB, r = 0) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.groundAnchorA.SetV(gaA);
    this.groundAnchorB.SetV(gaB);
    this.localAnchorA = this.bodyA.GetLocalPoint(anchorA);
    this.localAnchorB = this.bodyB.GetLocalPoint(anchorB);
    const d1X = anchorA.x - gaA.x;
    const d1Y = anchorA.y - gaA.y;
    this.lengthA = Math.sqrt(d1X * d1X + d1Y * d1Y);
    const d2X = anchorB.x - gaB.x;
    const d2Y = anchorB.y - gaB.y;
    this.lengthB = Math.sqrt(d2X * d2X + d2Y * d2Y);
    this.ratio = r;
    const C = this.lengthA + this.ratio * this.lengthB;
    this.maxLengthA = C - this.ratio * b2PulleyJoint.b2_minPulleyLength;
    this.maxLengthB = (C - b2PulleyJoint.b2_minPulleyLength) / this.ratio;
  }
};

const b2RevoluteJoint = class extends b2Joint {
  static tImpulse = new b2Vec2();
  constructor(def) {
    super(def);
    this.K = new b2Mat22();
    this.impulse3 = new b2Vec3();
    this.impulse2 = new b2Vec2();
    this.reduced = new b2Vec2();
    this.m_localAnchor1 = new b2Vec2();
    this.m_localAnchor2 = new b2Vec2();
    this.m_impulse = new b2Vec3();
    this.m_mass = new b2Mat33();
    this.m_localAnchor1.SetV(def.localAnchorA);
    this.m_localAnchor2.SetV(def.localAnchorB);
    this.m_referenceAngle = def.referenceAngle;
    this.m_impulse.SetZero();
    this.m_motorImpulse = 0.0;
    this.m_lowerAngle = def.lowerAngle;
    this.m_upperAngle = def.upperAngle;
    this.m_maxMotorTorque = def.maxMotorTorque;
    this.m_motorSpeed = def.motorSpeed;
    this.m_enableLimit = def.enableLimit;
    this.m_enableMotor = def.enableMotor;
    this.m_limitState = b2Joint.e_inactiveLimit;
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchor1);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchor2);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse.x, inv_dt * this.m_impulse.y);
  }
  GetReactionTorque(inv_dt = 0) {
    return inv_dt * this.m_impulse.z;
  }
  GetJointAngle() {
    return this.m_bodyB.m_sweep.a - this.m_bodyA.m_sweep.a - this.m_referenceAngle;
  }
  GetJointSpeed() {
    return this.m_bodyB.m_angularVelocity - this.m_bodyA.m_angularVelocity;
  }
  IsLimitEnabled() {
    return this.m_enableLimit;
  }
  EnableLimit(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableLimit = flag;
  }
  GetLowerLimit() {
    return this.m_lowerAngle;
  }
  GetUpperLimit() {
    return this.m_upperAngle;
  }
  SetLimits(lower = 0, upper = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_lowerAngle = lower;
    this.m_upperAngle = upper;
  }
  IsMotorEnabled() {
    return this.m_enableMotor;
  }
  EnableMotor(flag) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_enableMotor = flag;
  }
  SetMotorSpeed(speed = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_motorSpeed = speed;
  }
  GetMotorSpeed() {
    return this.m_motorSpeed;
  }
  SetMaxMotorTorque(torque = 0) {
    this.m_bodyA.SetAwake(true);
    this.m_bodyB.SetAwake(true);
    this.m_maxMotorTorque = torque;
  }
  GetMotorTorque() {
    return this.m_motorImpulse;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
    let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
    r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
    r1X = tX;
    tMat = bB.m_xf.R;
    let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
    let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
    r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
    r2X = tX;
    const m1 = bA.m_invMass;
    const m2 = bB.m_invMass;
    const i1 = bA.m_invI;
    const i2 = bB.m_invI;
    this.m_mass.col1.x = m1 + m2 + r1Y * r1Y * i1 + r2Y * r2Y * i2;
    this.m_mass.col2.x = (-r1Y * r1X * i1) - r2Y * r2X * i2;
    this.m_mass.col3.x = (-r1Y * i1) - r2Y * i2;
    this.m_mass.col1.y = this.m_mass.col2.x;
    this.m_mass.col2.y = m1 + m2 + r1X * r1X * i1 + r2X * r2X * i2;
    this.m_mass.col3.y = r1X * i1 + r2X * i2;
    this.m_mass.col1.z = this.m_mass.col3.x;
    this.m_mass.col2.z = this.m_mass.col3.y;
    this.m_mass.col3.z = i1 + i2;
    this.m_motorMass = i1 + i2;
    if (this.m_motorMass > 0.0) {
      this.m_motorMass = 1.0 / this.m_motorMass;
    }
    if (this.m_enableMotor == false) {
      this.m_motorImpulse = 0.0;
    }
    if (this.m_enableLimit) {
      const jointAngle = bB.m_sweep.a - bA.m_sweep.a - this.m_referenceAngle;
      if (b2Math.Abs(this.m_upperAngle - this.m_lowerAngle) < 2.0 * b2Settings.b2_angularSlop) {
        this.m_limitState = b2Joint.e_equalLimits;
      } else if (jointAngle <= this.m_lowerAngle) {
        if (this.m_limitState != b2Joint.e_atLowerLimit) {
          this.m_impulse.z = 0.0;
        }
        this.m_limitState = b2Joint.e_atLowerLimit;
      } else if (jointAngle >= this.m_upperAngle) {
        if (this.m_limitState != b2Joint.e_atUpperLimit) {
          this.m_impulse.z = 0.0;
        }
        this.m_limitState = b2Joint.e_atUpperLimit;
      } else {
        this.m_limitState = b2Joint.e_inactiveLimit;
        this.m_impulse.z = 0.0;
      }
    } else {
      this.m_limitState = b2Joint.e_inactiveLimit;
    }
    if (step.warmStarting) {
      this.m_impulse.x *= step.dtRatio;
      this.m_impulse.y *= step.dtRatio;
      this.m_impulse.z *= step.dtRatio;
      this.m_motorImpulse *= step.dtRatio;
      const PX = this.m_impulse.x;
      const PY = this.m_impulse.y;
      bA.m_linearVelocity.x -= m1 * PX;
      bA.m_linearVelocity.y -= m1 * PY;
      bA.m_angularVelocity -= i1 * ((r1X * PY - r1Y * PX) + this.m_motorImpulse + this.m_impulse.z);
      bB.m_linearVelocity.x += m2 * PX;
      bB.m_linearVelocity.y += m2 * PY;
      bB.m_angularVelocity += i2 * ((r2X * PY - r2Y * PX) + this.m_motorImpulse + this.m_impulse.z);
    } else {
      this.m_impulse.SetZero();
      this.m_motorImpulse = 0.0;
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat;
    let tX;
    let newImpulse;
    let r1X;
    let r1Y;
    let r2X;
    let r2Y;
    const v1 = bA.m_linearVelocity;
    let w1 = bA.m_angularVelocity;
    const v2 = bB.m_linearVelocity;
    let w2 = bB.m_angularVelocity;
    const m1 = bA.m_invMass;
    const m2 = bB.m_invMass;
    const i1 = bA.m_invI;
    const i2 = bB.m_invI;
    if (this.m_enableMotor && this.m_limitState != b2Joint.e_equalLimits) {
      const Cdot = w2 - w1 - this.m_motorSpeed;
      let impulse = this.m_motorMass * -Cdot;
      const oldImpulse = this.m_motorImpulse;
      const maxImpulse = step.dt * this.m_maxMotorTorque;
      this.m_motorImpulse = b2Math.Clamp(this.m_motorImpulse + impulse, -maxImpulse, maxImpulse);
      impulse = this.m_motorImpulse - oldImpulse;
      w1 -= i1 * impulse;
      w2 += i2 * impulse;
    }
    if (this.m_enableLimit && this.m_limitState != b2Joint.e_inactiveLimit) {
      tMat = bA.m_xf.R;
      r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
      r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
      r1X = tX;
      tMat = bB.m_xf.R;
      r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
      r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
      r2X = tX;
      const Cdot1X = v2.x + (-w2 * r2Y) - v1.x - (-w1 * r1Y);
      const Cdot1Y = v2.y + (w2 * r2X) - v1.y - (w1 * r1X);
      const Cdot2 = w2 - w1;
      this.m_mass.Solve33(this.impulse3, -Cdot1X, -Cdot1Y, -Cdot2);
      if (this.m_limitState == b2Joint.e_equalLimits) {
        this.m_impulse.Add(this.impulse3);
      } else if (this.m_limitState == b2Joint.e_atLowerLimit) {
        newImpulse = this.m_impulse.z + this.impulse3.z;
        if (newImpulse < 0.0) {
          this.m_mass.Solve22(this.reduced, -Cdot1X, -Cdot1Y);
          this.impulse3.x = this.reduced.x;
          this.impulse3.y = this.reduced.y;
          this.impulse3.z = -this.m_impulse.z;
          this.m_impulse.x += this.reduced.x;
          this.m_impulse.y += this.reduced.y;
          this.m_impulse.z = 0.0;
        }
      } else if (this.m_limitState == b2Joint.e_atUpperLimit) {
        newImpulse = this.m_impulse.z + this.impulse3.z;
        if (newImpulse > 0.0) {
          this.m_mass.Solve22(this.reduced, -Cdot1X, -Cdot1Y);
          this.impulse3.x = this.reduced.x;
          this.impulse3.y = this.reduced.y;
          this.impulse3.z = -this.m_impulse.z;
          this.m_impulse.x += this.reduced.x;
          this.m_impulse.y += this.reduced.y;
          this.m_impulse.z = 0.0;
        }
      }
      v1.x -= m1 * this.impulse3.x;
      v1.y -= m1 * this.impulse3.y;
      w1 -= i1 * (r1X * this.impulse3.y - r1Y * this.impulse3.x + this.impulse3.z);
      v2.x += m2 * this.impulse3.x;
      v2.y += m2 * this.impulse3.y;
      w2 += i2 * (r2X * this.impulse3.y - r2Y * this.impulse3.x + this.impulse3.z);
    } else {
      tMat = bA.m_xf.R;
      r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
      r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
      r1X = tX;
      tMat = bB.m_xf.R;
      r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
      r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
      r2X = tX;
      const CdotX = v2.x + (-w2 * r2Y) - v1.x - (-w1 * r1Y);
      const CdotY = v2.y + (w2 * r2X) - v1.y - (w1 * r1X);
      this.m_mass.Solve22(this.impulse2, -CdotX, -CdotY);
      this.m_impulse.x += this.impulse2.x;
      this.m_impulse.y += this.impulse2.y;
      v1.x -= m1 * this.impulse2.x;
      v1.y -= m1 * this.impulse2.y;
      w1 -= i1 * (r1X * this.impulse2.y - r1Y * this.impulse2.x);
      v2.x += m2 * this.impulse2.x;
      v2.y += m2 * this.impulse2.y;
      w2 += i2 * (r2X * this.impulse2.y - r2Y * this.impulse2.x);
    }
    bA.m_angularVelocity = w1;
    bB.m_angularVelocity = w2;
  }
  SolvePositionConstraints(baumgarte) {
    let C;
    let tMat;
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let angularError = 0.0;
    let positionError;
    let tX;
    let impulseX;
    let impulseY;
    if (this.m_enableLimit && this.m_limitState != b2Joint.e_inactiveLimit) {
      const angle = bB.m_sweep.a - bA.m_sweep.a - this.m_referenceAngle;
      let limitImpulse = 0.0;
      if (this.m_limitState == b2Joint.e_equalLimits) {
        C = b2Math.Clamp(angle - this.m_lowerAngle, -b2Settings.b2_maxAngularCorrection, b2Settings.b2_maxAngularCorrection);
        limitImpulse = -this.m_motorMass * C;
        angularError = b2Math.Abs(C);
      } else if (this.m_limitState == b2Joint.e_atLowerLimit) {
        C = angle - this.m_lowerAngle;
        angularError = -C;
        C = b2Math.Clamp(C + b2Settings.b2_angularSlop, -b2Settings.b2_maxAngularCorrection, 0.0);
        limitImpulse = -this.m_motorMass * C;
      } else if (this.m_limitState == b2Joint.e_atUpperLimit) {
        C = angle - this.m_upperAngle;
        angularError = C;
        C = b2Math.Clamp(C - b2Settings.b2_angularSlop, 0.0, b2Settings.b2_maxAngularCorrection);
        limitImpulse = -this.m_motorMass * C;
      }
      bA.m_sweep.a -= bA.m_invI * limitImpulse;
      bB.m_sweep.a += bB.m_invI * limitImpulse;
      bA.SynchronizeTransform();
      bB.SynchronizeTransform();
    }
    {
      tMat = bA.m_xf.R;
      let r1X = this.m_localAnchor1.x - bA.m_sweep.localCenter.x;
      let r1Y = this.m_localAnchor1.y - bA.m_sweep.localCenter.y;
      tX = tMat.col1.x * r1X + tMat.col2.x * r1Y;
      r1Y = tMat.col1.y * r1X + tMat.col2.y * r1Y;
      r1X = tX;
      tMat = bB.m_xf.R;
      let r2X = this.m_localAnchor2.x - bB.m_sweep.localCenter.x;
      let r2Y = this.m_localAnchor2.y - bB.m_sweep.localCenter.y;
      tX = tMat.col1.x * r2X + tMat.col2.x * r2Y;
      r2Y = tMat.col1.y * r2X + tMat.col2.y * r2Y;
      r2X = tX;
      let CX = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
      let CY = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
      const CLengthSquared = CX * CX + CY * CY;
      const CLength = Math.sqrt(CLengthSquared);
      positionError = CLength;
      const invMass1 = bA.m_invMass;
      const invMass2 = bB.m_invMass;
      const invI1 = bA.m_invI;
      const invI2 = bB.m_invI;
      const k_allowedStretch = 10.0 * b2Settings.b2_linearSlop;
      if (CLengthSquared > k_allowedStretch * k_allowedStretch) {
        const k = invMass1 + invMass2;
        const m = 1.0 / k;
        impulseX = m * -CX;
        impulseY = m * -CY;
        const k_beta = 0.5;
        bA.m_sweep.c.x -= k_beta * invMass1 * impulseX;
        bA.m_sweep.c.y -= k_beta * invMass1 * impulseY;
        bB.m_sweep.c.x += k_beta * invMass2 * impulseX;
        bB.m_sweep.c.y += k_beta * invMass2 * impulseY;
        CX = bB.m_sweep.c.x + r2X - bA.m_sweep.c.x - r1X;
        CY = bB.m_sweep.c.y + r2Y - bA.m_sweep.c.y - r1Y;
      }
      this.K.col1.x = invMass1 + invMass2 + invI1 * r1Y * r1Y + invI2 * r2Y * r2Y;
      this.K.col2.x = 0.0 + (-invI1 * r1X * r1Y) + (-invI2 * r2X * r2Y);
      this.K.col1.y = 0.0 + (-invI1 * r1X * r1Y) + (-invI2 * r2X * r2Y);
      this.K.col2.y = invMass1 + invMass2 + invI1 * r1X * r1X + invI2 * r2X * r2X;
      this.K.Solve(b2RevoluteJoint.tImpulse, -CX, -CY);
      impulseX = b2RevoluteJoint.tImpulse.x;
      impulseY = b2RevoluteJoint.tImpulse.y;
      bA.m_sweep.c.x -= bA.m_invMass * impulseX;
      bA.m_sweep.c.y -= bA.m_invMass * impulseY;
      bA.m_sweep.a -= bA.m_invI * (r1X * impulseY - r1Y * impulseX);
      bB.m_sweep.c.x += bB.m_invMass * impulseX;
      bB.m_sweep.c.y += bB.m_invMass * impulseY;
      bB.m_sweep.a += bB.m_invI * (r2X * impulseY - r2Y * impulseX);
      bA.SynchronizeTransform();
      bB.SynchronizeTransform();
    }
    return positionError <= b2Settings.b2_linearSlop && angularError <= b2Settings.b2_angularSlop;
  }
};

const b2RevoluteJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.type = b2Joint.e_revoluteJoint;
    this.localAnchorA.Set(0.0, 0.0);
    this.localAnchorB.Set(0.0, 0.0);
    this.referenceAngle = 0.0;
    this.lowerAngle = 0.0;
    this.upperAngle = 0.0;
    this.maxMotorTorque = 0.0;
    this.motorSpeed = 0.0;
    this.enableLimit = false;
    this.enableMotor = false;
  }
  Initialize(bA, bB, anchor) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA = this.bodyA.GetLocalPoint(anchor);
    this.localAnchorB = this.bodyB.GetLocalPoint(anchor);
    this.referenceAngle = this.bodyB.GetAngle() - this.bodyA.GetAngle();
  }
};

const b2WeldJoint = class extends b2Joint {
  constructor(def) {
    super(def);
    this.m_localAnchorA = new b2Vec2();
    this.m_localAnchorB = new b2Vec2();
    this.m_impulse = new b2Vec3();
    this.m_mass = new b2Mat33();
    this.m_tempImpulse = new b2Vec3();
    this.m_localAnchorA.SetV(def.localAnchorA);
    this.m_localAnchorB.SetV(def.localAnchorB);
    this.m_referenceAngle = def.referenceAngle;
    this.m_impulse.SetZero();
    this.m_mass = new b2Mat33();
  }
  GetAnchorA() {
    return this.m_bodyA.GetWorldPoint(this.m_localAnchorA);
  }
  GetAnchorB() {
    return this.m_bodyB.GetWorldPoint(this.m_localAnchorB);
  }
  GetReactionForce(inv_dt = 0) {
    return new b2Vec2(inv_dt * this.m_impulse.x, inv_dt * this.m_impulse.y);
  }
  GetReactionTorque(inv_dt = 0) {
    return inv_dt * this.m_impulse.z;
  }
  InitVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let rAX = this.m_localAnchorA.x - bA.m_sweep.localCenter.x;
    let rAY = this.m_localAnchorA.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rAX + tMat.col2.x * rAY;
    rAY = tMat.col1.y * rAX + tMat.col2.y * rAY;
    rAX = tX;
    tMat = bB.m_xf.R;
    let rBX = this.m_localAnchorB.x - bB.m_sweep.localCenter.x;
    let rBY = this.m_localAnchorB.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * rBX + tMat.col2.x * rBY;
    rBY = tMat.col1.y * rBX + tMat.col2.y * rBY;
    rBX = tX;
    const mA = bA.m_invMass;
    const mB = bB.m_invMass;
    const iA = bA.m_invI;
    const iB = bB.m_invI;
    this.m_mass.col1.x = mA + mB + rAY * rAY * iA + rBY * rBY * iB;
    this.m_mass.col2.x = (-rAY * rAX * iA) - rBY * rBX * iB;
    this.m_mass.col3.x = (-rAY * iA) - rBY * iB;
    this.m_mass.col1.y = this.m_mass.col2.x;
    this.m_mass.col2.y = mA + mB + rAX * rAX * iA + rBX * rBX * iB;
    this.m_mass.col3.y = rAX * iA + rBX * iB;
    this.m_mass.col1.z = this.m_mass.col3.x;
    this.m_mass.col2.z = this.m_mass.col3.y;
    this.m_mass.col3.z = iA + iB;
    if (step.warmStarting) {
      this.m_impulse.x *= step.dtRatio;
      this.m_impulse.y *= step.dtRatio;
      this.m_impulse.z *= step.dtRatio;
      bA.m_linearVelocity.x -= mA * this.m_impulse.x;
      bA.m_linearVelocity.y -= mA * this.m_impulse.y;
      bA.m_angularVelocity -= iA * (rAX * this.m_impulse.y - rAY * this.m_impulse.x + this.m_impulse.z);
      bB.m_linearVelocity.x += mB * this.m_impulse.x;
      bB.m_linearVelocity.y += mB * this.m_impulse.y;
      bB.m_angularVelocity += iB * (rBX * this.m_impulse.y - rBY * this.m_impulse.x + this.m_impulse.z);
    } else {
      this.m_impulse.SetZero();
    }
  }
  SolveVelocityConstraints(step) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    const vA = bA.m_linearVelocity;
    let wA = bA.m_angularVelocity;
    const vB = bB.m_linearVelocity;
    let wB = bB.m_angularVelocity;
    const mA = bA.m_invMass;
    const mB = bB.m_invMass;
    const iA = bA.m_invI;
    const iB = bB.m_invI;
    let tMat = bA.m_xf.R;
    let rAX = this.m_localAnchorA.x - bA.m_sweep.localCenter.x;
    let rAY = this.m_localAnchorA.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rAX + tMat.col2.x * rAY;
    rAY = tMat.col1.y * rAX + tMat.col2.y * rAY;
    rAX = tX;
    tMat = bB.m_xf.R;
    let rBX = this.m_localAnchorB.x - bB.m_sweep.localCenter.x;
    let rBY = this.m_localAnchorB.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * rBX + tMat.col2.x * rBY;
    rBY = tMat.col1.y * rBX + tMat.col2.y * rBY;
    rBX = tX;
    const Cdot1X = vB.x - wB * rBY - vA.x + wA * rAY;
    const Cdot1Y = vB.y + wB * rBX - vA.y - wA * rAX;
    const Cdot2 = wB - wA;
    const impulse = this.m_tempImpulse;
    this.m_mass.Solve33(impulse, -Cdot1X, -Cdot1Y, -Cdot2);
    this.m_impulse.Add(impulse);
    vA.x -= mA * impulse.x;
    vA.y -= mA * impulse.y;
    wA -= iA * (rAX * impulse.y - rAY * impulse.x + impulse.z);
    vB.x += mB * impulse.x;
    vB.y += mB * impulse.y;
    wB += iB * (rBX * impulse.y - rBY * impulse.x + impulse.z);
    bA.m_angularVelocity = wA;
    bB.m_angularVelocity = wB;
  }
  SolvePositionConstraints(baumgarte) {
    const bA = this.m_bodyA;
    const bB = this.m_bodyB;
    let tMat = bA.m_xf.R;
    let rAX = this.m_localAnchorA.x - bA.m_sweep.localCenter.x;
    let rAY = this.m_localAnchorA.y - bA.m_sweep.localCenter.y;
    let tX = tMat.col1.x * rAX + tMat.col2.x * rAY;
    rAY = tMat.col1.y * rAX + tMat.col2.y * rAY;
    rAX = tX;
    tMat = bB.m_xf.R;
    let rBX = this.m_localAnchorB.x - bB.m_sweep.localCenter.x;
    let rBY = this.m_localAnchorB.y - bB.m_sweep.localCenter.y;
    tX = tMat.col1.x * rBX + tMat.col2.x * rBY;
    rBY = tMat.col1.y * rBX + tMat.col2.y * rBY;
    rBX = tX;
    const mA = bA.m_invMass;
    const mB = bB.m_invMass;
    let iA = bA.m_invI;
    let iB = bB.m_invI;
    const C1X = bB.m_sweep.c.x + rBX - bA.m_sweep.c.x - rAX;
    const C1Y = bB.m_sweep.c.y + rBY - bA.m_sweep.c.y - rAY;
    const C2 = bB.m_sweep.a - bA.m_sweep.a - this.m_referenceAngle;
    const k_allowedStretch = 10.0 * b2Settings.b2_linearSlop;
    const positionError = Math.sqrt(C1X * C1X + C1Y * C1Y);
    const angularError = b2Math.Abs(C2);
    if (positionError > k_allowedStretch) {
      iA *= 1.0;
      iB *= 1.0;
    }
    this.m_mass.col1.x = mA + mB + rAY * rAY * iA + rBY * rBY * iB;
    this.m_mass.col2.x = (-rAY * rAX * iA) - rBY * rBX * iB;
    this.m_mass.col3.x = (-rAY * iA) - rBY * iB;
    this.m_mass.col1.y = this.m_mass.col2.x;
    this.m_mass.col2.y = mA + mB + rAX * rAX * iA + rBX * rBX * iB;
    this.m_mass.col3.y = rAX * iA + rBX * iB;
    this.m_mass.col1.z = this.m_mass.col3.x;
    this.m_mass.col2.z = this.m_mass.col3.y;
    this.m_mass.col3.z = iA + iB;
    const impulse = this.m_tempImpulse;
    this.m_mass.Solve33(impulse, -C1X, -C1Y, -C2);
    bA.m_sweep.c.x -= mA * impulse.x;
    bA.m_sweep.c.y -= mA * impulse.y;
    bA.m_sweep.a -= iA * (rAX * impulse.y - rAY * impulse.x + impulse.z);
    bB.m_sweep.c.x += mB * impulse.x;
    bB.m_sweep.c.y += mB * impulse.y;
    bB.m_sweep.a += iB * (rBX * impulse.y - rBY * impulse.x + impulse.z);
    bA.SynchronizeTransform();
    bB.SynchronizeTransform();
    return positionError <= b2Settings.b2_linearSlop && angularError <= b2Settings.b2_angularSlop;
  }
};

const b2WeldJointDef = class extends b2JointDef {
  constructor() {
    super();
    this.localAnchorA = new b2Vec2();
    this.localAnchorB = new b2Vec2();
    this.type = b2Joint.e_weldJoint;
    this.referenceAngle = 0.0;
  }
  Initialize(bA, bB, anchor) {
    this.bodyA = bA;
    this.bodyB = bB;
    this.localAnchorA.SetV(this.bodyA.GetLocalPoint(anchor));
    this.localAnchorB.SetV(this.bodyB.GetLocalPoint(anchor));
    this.referenceAngle = this.bodyB.GetAngle() - this.bodyA.GetAngle();
  }
};

export const Box2D = {
  Collision: {
    b2AABB,
    b2Bound,
    b2BoundValues,
    b2Collision,
    b2ContactID,
    b2ContactPoint,
    b2Distance,
    b2DistanceInput,
    b2DistanceOutput,
    b2DistanceProxy,
    b2DynamicTree,
    b2DynamicTreeBroadPhase,
    b2DynamicTreeNode,
    b2DynamicTreePair,
    b2Manifold,
    b2ManifoldPoint,
    b2Point,
    b2RayCastInput,
    b2RayCastOutput,
    b2Segment,
    b2SeparationFunction,
    b2Simplex,
    b2SimplexCache,
    b2SimplexVertex,
    b2TimeOfImpact,
    b2TOIInput,
    b2WorldManifold,
    ClipVertex,
    Features,
    Shapes: {
      b2CircleShape,
      b2EdgeChainDef,
      b2EdgeShape,
      b2MassData,
      b2PolygonShape,
      b2Shape
    }
  },
  Common: {
    b2Color,
    b2Settings,
    Math: {
      b2Mat22,
      b2Mat33,
      b2Math,
      b2Sweep,
      b2Transform,
      b2Vec2,
      b2Vec3
    }
  },
  Dynamics: {
    b2Body,
    b2BodyDef,
    b2ContactFilter,
    b2ContactImpulse,
    b2ContactListener,
    b2ContactManager,
    b2DebugDraw,
    b2DestructionListener,
    b2FilterData,
    b2Fixture,
    b2FixtureDef,
    b2Island,
    b2TimeStep,
    b2World,
    Contacts: {
      b2CircleContact,
      b2Contact,
      b2ContactConstraint,
      b2ContactConstraintPoint,
      b2ContactEdge,
      b2ContactFactory,
      b2ContactRegister,
      b2ContactResult,
      b2ContactSolver,
      b2EdgeAndCircleContact,
      b2NullContact,
      b2PolyAndCircleContact,
      b2PolyAndEdgeContact,
      b2PolygonContact,
      b2PositionSolverManifold
    },
    Controllers: {
      b2BuoyancyController,
      b2ConstantAccelController,
      b2ConstantForceController,
      b2Controller,
      b2ControllerEdge,
      b2GravityController,
      b2TensorDampingController
    },
    Joints: {
      b2DistanceJoint,
      b2DistanceJointDef,
      b2FrictionJoint,
      b2FrictionJointDef,
      b2GearJoint,
      b2GearJointDef,
      b2Jacobian,
      b2Joint,
      b2JointDef,
      b2JointEdge,
      b2LineJoint,
      b2LineJointDef,
      b2MouseJoint,
      b2MouseJointDef,
      b2PrismaticJoint,
      b2PrismaticJointDef,
      b2PulleyJoint,
      b2PulleyJointDef,
      b2RevoluteJoint,
      b2RevoluteJointDef,
      b2WeldJoint,
      b2WeldJointDef
    }
  }
};
