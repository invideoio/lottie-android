package com.airbnb.lottie

import android.graphics.Path
import android.graphics.Point
import android.graphics.PointF
import android.os.Build
import androidx.annotation.RequiresApi
import com.soywiz.korma.geom.bezier.Bezier
import com.soywiz.korma.geom.bezier.SegmentEmitter
import com.soywiz.korma.geom.vector.VectorPath
import com.soywiz.korma.geom.vector.isEmpty
import com.soywiz.korma.geom.vector.lineTo
import kotlin.math.*

const val kMaxTValue = 0x3FFFFFFF
const val CHEAP_DIST_LIMIT = 0.5F

abstract class IVPathMeasure(path: VectorPath = VectorPath(), forceClosed : Boolean = false,
                             resScale: Float) {

    var path: VectorPath
    var forceClosed: Boolean
    var isClosed: Boolean = false
    var segments = mutableListOf<Segment>()
    var pts = mutableListOf<PointF>()
    var verbs = mutableListOf<Int>()
    var firstPtIndex : Int = -1
    private var length = -1F
    lateinit var iter: Iter
    var tolerance : Float = CHEAP_DIST_LIMIT * resScale

    init {
        this.path = path
        this.forceClosed = forceClosed
    }

    enum class Verb {
        kMove_Verb,
        kLine_Verb,
        kQuad_Verb,
        kConic_Verb,
        kCubic_Verb,
        kClose_Verb,
        kDone_Verb
    }

    class Iter {

    }

    fun Iter.setPath(path: VectorPath, forceClosed: Boolean) {
        pts.addAll(path.getPoints() as MutableList<PointF>)
        verbs.addAll(path.commands)
        // fVerbStop not declared or set since it won't be necessary for now
        // fConicWeights not declared or set since Conics are not being handled yet

        // some variables of SkPath do not exist on VectorPath class, are they gonna be
        // necessary ???

    }

    fun Iter.next(pts: Array<PointF>, doConsumeDegenerates: Boolean = true,
             exact: Boolean = false) : Verb {
        if(doConsumeDegenerates) {
            this.consumeDegenerateSegments(exact)
        }
        return this.doNext(pts)
    }

    fun Iter.consumeDegenerateSegments(exact: Boolean) {
        // MUST Implement
    }

    fun Iter.doNext(pts: Array<PointF>) : Verb {
        return Verb.kLine_Verb // MUST Implement properly
    }

    companion object {
        fun tValue2Scalar(t: Int): Float {
            return t.toFloat() / kMaxTValue
        }
    }

    data class Segment(var distance: Float, var ptIndex: Int,
                       var tValue: Int = 30, var type: SegType = SegType.kQuad_SegType) {
        fun getT() : Float {
            return tValue2Scalar(this.tValue)
        }
    }

    data class QuadCoeff(var A: PointF, var B: PointF, var C: PointF) {

        fun eval(tt: PointF) : PointF {
            return PointF((this.A.x * tt.x + this.B.x) * tt.x + this.C.x,
                          (this.A.y * tt.y + this.B.y) * tt.y + this.C.y)
        }
    }

    fun QuadCoeff(src: Array<PointF>) : QuadCoeff {
        return QuadCoeff(src[0], src[1], src[2])
    }

    data class CubicCoeff(var A: PointF, var B: PointF, var C: PointF, var D: PointF) {
        fun eval(t: Float) : PointF {
            return PointF(((A.x * t + B.x) * t + C.x) * t + D.x,
                          ((A.y * t + B.y) * t + C.y) * t + D.y)
        }
    }

    fun CubicCoeff(src: Array<PointF>): CubicCoeff {
        return CubicCoeff(src[0], src[1], src[2], src[3])
    }

    data class Vec2(var x: Float, var y: Float) {
        fun normalize(x: Float, y: Float) {
            val length = PointF.length(x, y)
            this.x = x / length
            this.y = y / length
        }

        fun set(pt: Vec2) {
            this.x = pt.x
            this.y = pt.y
        }
    }

    private fun Vec2(point: PointF): IVPathMeasure.Vec2 {
        return Vec2(point.x, point.y)
    }

    enum class SegType {
        kLine_SegType,
        kQuad_SegType,
        kCubic_SegType,
        kConic_SegType,
    }

    fun lerp(a: Float, b: Float, t: Float) : Float {
        if(t <= 0.0f) {
            return a
        } else if(t >= 1.0f) {
            return b
        } else {
            return a + (b - a) * t
        }
    }

    fun lerp(a: PointF, b: PointF, t: PointF) : PointF {

        var x = 0F
        if(t.x <= 0.0f) {
            x = a.x
        } else if(t.x >= 1.0f) {
            x = b.x
        } else {
            x = a.x + (b.x - a.x) * t.x
        }

        var y = 0F
        if(t.y <= 0.0f) {
            y = a.y
        } else if(t.y >= 1.0f) {
            x = b.y
        } else {
            x = a.y + (b.y - a.y) * t.y
        }

        return PointF(x, y)
    }

    private fun tspan_big_enough(tspan: Int): Boolean {
        return (tspan shr 10) > 0
    }

    fun quad_too_curvy(pts: Array<PointF>) : Boolean {
        // diff = (a/4 + b/2 + c/4) - (a/2 + c/2)
        // diff = -a/4 + b/2 - c/4
        val dx = 0.5f * pts[1].x -
                 0.25f * (pts[0].x + pts[2].x)
        val dy = 0.5f * pts[1].y -
                0.25f * (pts[0].y + pts[2].y)

        val dist = max(dx, dy)
        return dist > tolerance
    }

//    fun conic_too_curvy(firstPt: PointF, midTPt: PointF, lastPt: PointF) : Boolean {
//        val midEnds : PointF = PointF(0.5f * (firstPt.x + lastPt.x),
//                                      0.5f * (firstPt.y + lastPt.y))
//        val dxy : PointF = PointF(midTPt.x - midEnds.x,
//                                  midTPt.y - midEnds.y)
//        val dist = max(abs(dxy.x), abs(dxy.y))
//        return dist > tolerance
//    }

    fun cheap_dist_exceeds_limits(pt: PointF, x: Float, y: Float) : Boolean {
        val dist = max(abs(x - pt.x), abs(y - pt.y))
        return dist > tolerance
    }

    fun valid_unit_divide(numerator: Float, denominator: Float) : Float {
        var numer = numerator
        var denom = denominator
        if(numer < 0F) {
            numer = -numer
            denom = -denom
        }

        if(denom == 0F || numer == 0F || numer >= denom) {
            return 0F
        }

        val r = numer / denom

        if(r == Float.NaN) {
            return 0F
        }

        if(r == 0F) {
            return 0F
        }

        return r
    }

    fun findQuadMaxCurvature(src : Array<PointF>) : Float {
        val Ax = src[1].x - src[0].x
        val Ay = src[1].y - src[0].y
        val Bx = src[0].x - src[1].x - src[1].x + src[2].x
        val By = src[0].y - src[1].y - src[1].y + src[2].y
        val t = valid_unit_divide(-(Ax * Bx + Ay * By), Bx * Bx + By * By);
        return t;
    }

    fun cubic_too_curvy(pts: Array<PointF>) : Boolean {
        return cheap_dist_exceeds_limits(pts[1],
                                         lerp(pts[0].x, pts[3].x, 1.0f / 3.0f),
                                         lerp(pts[0].y, pts[3].y, 1.0f / 3.0f)) ||
               cheap_dist_exceeds_limits(pts[2],
                                         lerp(pts[0].x, pts[3].x, 2.0f / 3.0f),
                                         lerp(pts[0].y, pts[3].y, 2.0f / 3.0f))
    }

    fun evalQuadAt(src: Array<PointF>, t: Float) : PointF {
        val tt = PointF(t, t)
        return QuadCoeff(src).eval(tt)
    }

    fun evalQuadAt(src: Array<PointF>, t: Float, pt: PointF, tangent: Vec2) {
        pt.set(evalQuadAt(src, t))
        tangent.set(evalQuadTangentAt(src, t))
    }

    private fun evalQuadTangentAt(src: Array<PointF>, t: Float): Vec2 {
        val p0 = src[0]
        val p1 = src[1]
        val p2 = src[2]

        val B = PointF(p1.x - p0.x,
                       p1.y - p1.y)
        val A = PointF(p2.x - p1.x - B.x,
                       p2.y - p1.y - B.y)

        return Vec2(A.x * t + B.x, A.y * t + B.y)
    }

    fun quad_folden_len(pts: Array<PointF>) : Float {
        val t = findQuadMaxCurvature(pts)
        val pt = evalQuadAt(pts, t)
        val a = PointF(pts[2].x - pt.x,
                       pts[2].y - pt.y)
        var result = a.length()
        if(0F != t) {
            val b = PointF(pts[0].x - pt.x,
                           pts[0].y - pt.y)
            result += b.length()
        }
        return result
    }

    fun compute_quad_len(pts: Array<PointF>) : Float {
        var a = PointF(0F,0F)
        var b = PointF(0F,0F)
        a.x = pts[0].x - 2F * pts[1].x + pts[2].x
        a.y = pts[0].y - 2F * pts[1].y + pts[2].y
        val A = 4F * (a.x * a.x + a.y * a.y)
        if(A == 0F) {
            a.x = pts[2].x - pts[0].x
            a.y = pts[2].y - pts[0].y
            return a.length()
        }
        b.x = 2F * (pts[1].x - pts[0].x)
        b.y = 2F * (pts[1].y - pts[0].y)
        val B = 4F * (a.x * b.x + a.y * b.y)
        val C =      (b.x * b.x + b.y * b.y)
        val Sabc = 2F * sqrt(A + B + C)
        val A_2 = sqrt(A)
        val A_32 = 2F * A * A_2
        val C_2 = 2F * sqrt(C)
        val BA = B / A_2
        if(BA + C_2 == 0F) {
            return quad_folden_len(pts)
        }
        val J = A_32 * Sabc + A_2 * B * (Sabc - C_2)
        val K = 4F * C * A - B * B
        val L = (2F * A_2 + BA + Sabc) / (BA + C_2)
        if(L <= 0F) {
            return quad_folden_len(pts)
        }
        val M = ln(L)
        return (J + K * M) / (4F * A_32)
    }

    private fun chopQuadAt(pts: Array<PointF>, dst: Array<PointF>, t: Float) {
        val p0 = pts[0]
        val p1 = pts[1]
        val p2 = pts[2]
        val tt = PointF(t, t)

        val p01 = lerp(p0, p1, tt)
        val p12 = lerp(p1, p2, tt)

        dst[0] = p0;
        dst[0] = p01;
        dst[0] = lerp(p01, p12, tt)
        dst[0] = p12;
        dst[0] = p2;
    }

    private fun chopQuadAtHalf(pts: Array<PointF>, dst: Array<PointF>) {
        chopQuadAt(pts, dst, 0.5f)
    }

    fun compute_quad_segs(pts: Array<PointF>, distance: Float,
                          mint: Int, maxt: Int, ptIndex: Int) : Float {
        var newDist = 0F
        if(tspan_big_enough(maxt - mint) && quad_too_curvy(pts)) {
            val tmp = Array(5) { PointF() }
            val halfT = (mint + maxt) shr 1

            chopQuadAtHalf(pts, tmp)
            newDist = distance

            val first_tmp = arrayOf(tmp[0], tmp[1], tmp[2])
            val second_tmp = arrayOf(tmp[2], tmp[3], tmp[4])
            newDist = compute_quad_segs(first_tmp, distance, mint, halfT, ptIndex)
            newDist = compute_quad_segs(second_tmp, newDist, halfT, maxt, ptIndex)
        } else {
            val d = PointF.length(pts[0].x - pts[2].x, pts[0].y - pts[2].y)
            val prevD = newDist
            newDist += d
            if(newDist > prevD) {
                segments.add(Segment(newDist, ptIndex, maxt, SegType.kQuad_SegType))
            }
        }
        return newDist
    }

//    fun compute_conic_segs(conic: Conic, distance: Float, mint: Int, minPt: PointF,
//                           maxt: Int, maxPt: PointF, ptIndex: Int) : Float {
//        val halfT = (mint + maxt) shr 1
//        val halfPt = conic.evalAt(tValue2Scalar(halfT))
//        var newDist = 0F
//        if(tspan_big_enough(maxt - mint) && conic_too_curvy(minPt, halfPt, maxPt)) {
//            newDist = compute_conic_segs(conic, distance, mint, minPt, halfT, halfPt, ptIndex)
//            newDist = compute_conic_segs(conic, newDist, halfT, halfPt, maxt, maxPt, ptIndex)
//        } else {
//            val d = PointF.length(minPt.x - maxPt.x, minPt.y - maxPt.y)
//            val prevD = distance
//            newDist = distance + d
//            if(newDist > prevD) {
//                segments.add(Segment(newDist, ptIndex, maxt, SegType.kQuad_SegType))
//            }
//        }
//        return newDist
//    }

    private fun chopCubicAt(src: Array<PointF>, dst: Array<PointF>, t: Float) {
        val p0 = src[0]
        val p1 = src[1]
        val p2 = src[2]
        val p3 = src[3]
        val tt = PointF(t, t)

        val ab = lerp(p0, p1, tt)
        val bc = lerp(p1, p2, tt)
        val cd = lerp(p2, p3, tt)
        val abc = lerp(ab, bc, tt)
        val bcd = lerp(bc, cd, tt)
        val abcd = lerp(abc, bcd, tt)

        dst[0] = src[0]
        dst[1] = ab
        dst[2] = abc
        dst[3] = abcd
        dst[4] = bcd
        dst[5] = cd
        dst[6] = src[3]
    }

    private fun chopCubicAtHalf(src: Array<PointF>, dst: Array<PointF>) {
        chopCubicAt(src, dst, 0.5f)
    }

    fun compute_cubic_segs(pts: Array<PointF>, distance: Float, mint: Int,
                           maxt: Int, ptIndex: Int) : Float {
        var newDist = 0F
        if(tspan_big_enough(maxt - mint) && cubic_too_curvy(pts)) {
            val tmp = Array(7) { PointF() }
            val halfT = (mint + maxt) shr 1
            chopCubicAtHalf(pts, tmp)

            val first_tmp = arrayOf(tmp[0], tmp[1], tmp[2], tmp[3])
            val second_tmp = arrayOf(tmp[3], tmp[4], tmp[5], tmp[6])
            newDist = compute_cubic_segs(first_tmp, distance, mint, halfT, ptIndex)
            newDist = compute_cubic_segs(second_tmp, newDist, halfT, maxt, ptIndex)
        } else {
            val d = PointF.length(pts[0].x - pts[3].x, pts[0].y - pts[3].y)
            val prevD = newDist
            newDist += d
            if(newDist > prevD) {
                segments.add(Segment(newDist, ptIndex, maxt, SegType.kCubic_SegType))
            }
        }
        return newDist
    }

    fun segTo(pts: Array<PointF>, segType: SegType, startT: Float, stopT: Float, dst: VectorPath) {
        if(startT == stopT) {
            val lastPt = dst.getPoints().last() as PointF
            dst.lineTo(lastPt.x, lastPt.y)
            return
        }

        var tmp0 = Array(7) { PointF() }
        var tmp1 = Array(7) { PointF() }

//        val conic = Conic(pts[0], pts[2], pts[3], pts[1].x)

        when(segType) {
            // Straight segments
            SegType.kLine_SegType ->
                if(stopT == 1F) {
                    dst.lineTo(pts[1].x.toDouble(),
                               pts[1].y.toDouble())
                } else {
                    val stopTT = PointF(stopT, stopT)
                    val pts01 = lerp(pts[0], pts[1], stopTT)
                    dst.lineTo(pts01.x.toDouble(),
                               pts01.y.toDouble())
                }

            // Quadratic segments
            SegType.kQuad_SegType ->
                if(startT == 0F) {
                    if(stopT == 1F) {
                        dst.quadTo(pts[1].x.toDouble(), pts[1].y.toDouble(),
                                   pts[2].x.toDouble(), pts[2].y.toDouble())
                    } else {
                        chopQuadAt(pts, tmp0, stopT)
                        dst.quadTo(tmp0[1].x.toDouble(), tmp0[1].y.toDouble(),
                                   tmp0[2].x.toDouble(), tmp0[2].y.toDouble())
                    }
                } else {
                    chopQuadAt(pts, tmp0, startT)
                    if(stopT == 1F) {
                        dst.quadTo(tmp0[3].x.toDouble(), tmp0[3].y.toDouble(),
                                   tmp0[4].x.toDouble(), tmp0[4].y.toDouble())
                    } else {
                        val second_tmp = arrayOf(tmp0[2], tmp0[3], tmp0[4])
                        chopQuadAt(second_tmp, tmp1, (stopT - startT) / (1F - startT))
                        dst.quadTo(tmp1[1].x.toDouble(), tmp1[1].y.toDouble(),
                                   tmp1[2].x.toDouble(), tmp1[2].y.toDouble())
                    }
                }

            // Conic segments
            SegType.kConic_SegType -> println("Conic SegType not implemented")
//                if(startT == 0F) {
//                    if(stopT == 1F) {
//                        dst.conicTo(conic.pts[1], conic.pts[2], conic.w)
//                    } else {
//                        val tmp : List<Conic>
//                        if(conic.chopAt(stopT, tmp1)) {
//                            dst.conicTo(tmp[0].pts[1], tmp[0].pts[2], tmp[0].w)
//                        }
//                    }
//                } else {
//                    if(stopT == 1F) {
//                        val tmp : List<Conic>
//                        if(conic.chopAt(startT, tmp)) {
//                            dst.conicTo(tmp[1].pts[1], tmp[1].pts[2], tmp[1].w)
//                        }
//                    } else {
//                        val tmp : List<Conic>
//                        conic.chopAt(startT, stopT, tmp)
//                        dst.conicTo(tmp[0].pts[1], tmp[0].pts[2], tmp[0].w)
//                    }
//                }

            // Cubic segments
            SegType.kCubic_SegType ->
                if(startT == 0F) {
                    if(stopT == 1F) {
                        dst.cubicTo(pts[1].x.toDouble(), pts[1].y.toDouble(),
                                    pts[2].x.toDouble(), pts[2].y.toDouble(),
                                    pts[3].x.toDouble(), pts[3].y.toDouble())
                    } else {
                        chopCubicAt(pts, tmp0, stopT)
                        dst.cubicTo(tmp0[1].x.toDouble(), tmp0[1].y.toDouble(),
                                    tmp0[2].x.toDouble(), tmp0[2].y.toDouble(),
                                    tmp0[3].x.toDouble(), tmp0[3].y.toDouble())
                    }
                } else {
                    chopCubicAt(pts, tmp0, startT)
                    if(stopT == 1F) {
                        dst.cubicTo(tmp0[4].x.toDouble(), tmp0[4].y.toDouble(),
                                    tmp0[5].x.toDouble(), tmp0[5].y.toDouble(),
                                    tmp0[6].x.toDouble(), tmp0[6].y.toDouble())
                    } else {
                        val second_tmp = arrayOf(tmp0[3], tmp0[4], tmp0[5], tmp0[6])
                        chopCubicAt(second_tmp, tmp1, (stopT - startT) / (1F - startT))
                        dst.cubicTo(tmp1[1].x.toDouble(), tmp1[1].y.toDouble(),
                                    tmp1[2].x.toDouble(), tmp1[2].y.toDouble(),
                                    tmp1[3].x.toDouble(), tmp1[3].y.toDouble())
                    }
                }

        }

    }

    fun compute_pos_tan(pts: Array<PointF>, segType: SegType, t: Float,
                        pos: PointF, tangent: Vec2) {
        when(segType) {
            SegType.kLine_SegType -> {
                val tt = PointF(t, t)
                pos.set(lerp(pts[0], pts[1], tt))
                tangent.normalize(pts[1].x - pts[0].x, pts[1].y - pts[0].y)
            }
            SegType.kQuad_SegType -> {
                evalQuadAt(pts, t, pos, tangent)
                tangent.normalize(tangent.x, tangent.y)
            }
            SegType.kConic_SegType -> println("Conic SegType not implemented")
//            {
//                Conic(pts[0], pts[2], pts[3], pts[1].x).evalAt(t, pos, tangent)
//                if (tangent.isNotEmpty()) {
//                    tangent[0].normalize(tangent[0].x, tangent[0].y)
//                }
//            }
            SegType.kCubic_SegType -> {
                evalCubicAt(pts, t, pos, tangent) // NOTE: Originally, this receives a nullptr,
                                                  // is there anywhere on code where this final
                                                  // value is passed ???
                tangent.normalize(tangent.x, tangent.y)
            }
        }
    }

    private fun evalCubicAt(src: Array<PointF>, t: Float,
                            loc: PointF, tangent: IVPathMeasure.Vec2) {
        loc.set(CubicCoeff(src).eval(t))
        if((t == 0F && src[0] == src[1]) || (t == 1F && src[2] == src[3])) {
            if(t == 0F) {
                val vec02 = Vec2(src[2].x - src[0].x,
                                 src[2].y - src[0].y)
                tangent.set(vec02)
            } else {
                val vec13 = Vec2(src[3].x - src[1].x,
                                 src[3].y - src[1].y)
                tangent.set(vec13)
            }
        } else {
            tangent.set(eval_cubic_derivative(src, t))
        }

        // NOTE: This function has a 4th input pointer called curvature, which is not being used
        // anywhere in the code...
    }

    private fun eval_cubic_derivative(src: Array<PointF>, t: Float): IVPathMeasure.Vec2 {
        val coeff = QuadCoeff(PointF(), PointF(), PointF())
        val p0 = src[0]
        val p1 = src[1]
        val p2 = src[2]
        val p3 = src[3]

        coeff.A = PointF(p3.x + 3F * (p1.x - p2.x) - p0.x,
                         p3.y + 3F * (p1.y - p2.y) - p0.y)
        coeff.B = PointF(2F * (p2.x - 2F * p1.x + p0.x),
                         2F * (p2.y - 2F * p1.y + p0.y))
        coeff.C = PointF(p1.x - p0.x,
                         p1.y - p0.y)

        val tt = PointF(t, t)
        return Vec2(coeff.eval(tt))
    }

    fun buildSegments() {
        val pts = Array(4){ PointF() }
        var ptIndex = this.firstPtIndex
        var distance = 0F
        var isClosed = this.forceClosed
        var firstMoveTo = ptIndex < 0
        val seg = Segment(0F, 0)

        segments.clear()
        var done = false
        do {
            when(iter.next(pts)) {
                Verb.kMove_Verb -> {
                    ptIndex += 1
                    this.pts.add(pts[0])
                    if(!firstMoveTo) {
                        done = true
                        break
                    }
                    firstMoveTo = false
                    break
                }
                Verb.kLine_Verb -> {
                    val d = PointF.length(pts[0].x - pts[1].x, pts[0].y - pts[1].y)
                    val prevD = distance
                    distance += d
                    if(distance > prevD) {
                        segments.add(Segment(distance, ptIndex, kMaxTValue, SegType.kLine_SegType))
                        this.pts.add(pts[1])
                        ptIndex++
                    }
                    break
                }
                Verb.kQuad_Verb -> {
                    val prevD = distance
                    // NOTE: There was and if(false) here originally, which will never trigger.
                    // Why as this here ???
                    distance = compute_quad_segs(pts, distance, 0, kMaxTValue, ptIndex)
                    if(distance > prevD) {
                        this.pts.add(pts[1])
                        this.pts.add(pts[2])
                        ptIndex += 2
                    }
                    break
                }
                Verb.kConic_Verb -> println("Conic SegType not implemented")
                Verb.kCubic_Verb -> {
                    val prevD = distance
                    distance = compute_cubic_segs(pts, distance, 0, kMaxTValue, ptIndex)
                    if(distance > prevD) {
                        this.pts.add(pts[1])
                        this.pts.add(pts[2])
                        this.pts.add(pts[3])
                        ptIndex += 3
                    }
                    break
                }
                Verb.kClose_Verb -> {
                    isClosed = true
                    break
                }
                Verb.kDone_Verb -> {
                    done = true
                    break
                }
            }
        } while(!done)

        this.length = distance
        this.isClosed = isClosed
        this.firstPtIndex = ptIndex
    }

    fun length() : Float {
        if(this.path.isEmpty()) {
            return 0F
        }
        if(this.length < 0F) {
            buildSegments()
        }
        if(length.isNaN()) {
            this.length = 0F
        }
        return this.length
    }

    fun nextContour() : Boolean {
        this.length = -1F
        return length() > 0F
    }

    fun setPath(path: VectorPath, forceClosed: Boolean) {
        this.path = path
        this.length = -1F   // signal we need to compute it
        this.forceClosed = forceClosed
        this.firstPtIndex = -1
        this.iter.setPath(path, forceClosed)
        this.segments.clear()
        this.pts.clear()
    }

}