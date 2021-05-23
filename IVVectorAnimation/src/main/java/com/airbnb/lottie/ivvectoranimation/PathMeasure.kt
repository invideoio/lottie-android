@file:Suppress("unused")

package com.airbnb.lottie.ivvectoranimation

import android.graphics.Matrix
import android.graphics.PointF
import android.graphics.RectF
import android.util.SparseArray
import com.soywiz.korma.geom.BoundsBuilder
import com.soywiz.korma.geom.Rectangle
import com.soywiz.korma.geom.vector.*
import java.util.*
import kotlin.math.*

const val kMaxTValue = 0x3FFFFFFF
const val CHEAP_DIST_LIMIT = 0.5F

class Path {

    lateinit var bounds: RectF
        private set

    var fillType: FillType = FillType.EVEN_ODD
    var path: VectorPath = VectorPath()

    fun moveTo(x: Float, y: Float) {
        path.moveTo(x, y)
    }

    fun quadTo(fl: Float, fl1: Float, fl2: Float, fl3: Float) {
        path.quadTo(fl, fl1, fl2, fl3)
    }

    fun cubicTo(fl: Float, fl1: Float, fl2: Float, fl3: Float, x: Float, y: Float) {
        path.cubicTo(fl, fl1, fl2, fl3, x, y)
    }

    fun lineTo(x: Float, y: Float) {
        path.lineTo(x, y)
    }

    fun reset() {
        path.clear()
    }

    fun addPath(tempPath: Path) {
        path.appendFrom(tempPath.path)
    }

    fun set(tempPath: Path) {
        path.setFrom(tempPath.path)
    }

    fun close() {
        path.close()
    }

    fun offset(x: Float, y: Float) {
        val matrix = com.soywiz.korma.geom.Matrix(
            1F, 0F,
            0F, 1F,
            x, y
        )
        path.applyTransform(matrix)
    }

    fun addPath(tempPath: Path, matrix: Matrix) {
        val aPath = Path()
        aPath.set(tempPath)
        aPath.transform(matrix)
        path.appendFrom(aPath.path)
    }

    fun computeBounds(rect: RectF, b: Boolean) {
        val rectangle = Rectangle()
        path.getBounds(rectangle, BoundsBuilder())
        rect.set(
            rectangle.left.toFloat(),
            rectangle.top.toFloat(),
            rectangle.right.toFloat(),
            rectangle.bottom.toFloat()
        )
        bounds = RectF(rect)
    }

    fun transform(matrix: Matrix) {
        val values = FloatArray(9)
        matrix.getValues(values)
        val korgeMatrix = com.soywiz.korma.geom.Matrix(
            values[0], values[1],
            values[3], values[4],
            values[2], values[5]
        )
        path.applyTransform(korgeMatrix)
    }

    fun arcTo(rect: RectF, i: Int, i1: Int, b: Boolean) {
        // TODO: Properly implement it later
    }

    /**
     * Set this path to the result of applying the Op to this path and the specified path.
     * The resulting path will be constructed from non-overlapping contours.
     * The curve order is reduced where possible so that cubics may be turned
     * into quadratics, and quadratics maybe turned into lines.
     *
     * @param path The second operand (for difference, the subtrahend)
     *
     * @return True if operation succeeded, false otherwise and this path remains unmodified.
     *
     * @see Op
     * @see #op(Path, Path, android.graphics.Path.Op)
     */
    open fun op(path: Path, op: Op) : Boolean {
        return op(this, path, op)
    }


    /**
     * Set this path to the result of applying the Op to the two specified paths.
     * The resulting path will be constructed from non-overlapping contours.
     * The curve order is reduced where possible so that cubics may be turned
     * into quadratics, and quadratics maybe turned into lines.
     *
     * @param path1 The first operand (for difference, the minuend)
     * @param path2 The second operand (for difference, the subtrahend)
     *
     * @return True if operation succeeded, false otherwise and this path remains unmodified.
     *
     * @see Op
     * @see #op(Path, android.graphics.Path.Op)
     */

    open fun op(path1: Path, path2: Path, op: Op): Boolean {
        // Op(*p1, *p2, op, r);
/*        if (nOp(path1.mNativePath, path2.mNativePath, op.ordinal(), this.mNativePath)) {
            isSimplePath = false
            rects = null
            return true
        }*/
        return false
    }


    enum class FillType {
        EVEN_ODD,
        INVERSE_EVEN_ODD,
        INVERSE_WINDING,
        WINDING
    }

    /**
     * The logical operations that can be performed when combining two paths.
     *
     * @see #op(Path, android.graphics.Path.Op)
     * @see #op(Path, Path, android.graphics.Path.Op)
     */
    enum class Op {
        /**
         * Subtract the second path from the first path.
         */
        DIFFERENCE,
        /**
         * Intersect the two paths.
         */
        INTERSECT,
        /**
         * Union (inclusive-or) the two paths.
         */
        UNION,
        /**
         * Exclusive-or the two paths.
         */
        XOR,
        /**
         * Subtract the first path from the second path.
         */
        REVERSE_DIFFERENCE
    }

// THIS FUNCTION IS NOT BEING USED IN PRODUCTION
//    fun computeBounds(bounds: RectF, b: Boolean) {
//
//    }

}

class PathMeasure(
    var path: VectorPath = VectorPath(),
    private var forceClosed: Boolean = false,
    private var resScale: Float = 1F
) {
    private var isClosed: Boolean = false
    private var segments = mutableListOf<Segment>()
    private var pts = mutableListOf<PointF>()
    private var firstPtIndex: Int = -1
    private var pathLength = -1F
    var tolerance : Float
    val arrayPool = SparseArray<Queue<Array<PointF>>>()
    init {
        tolerance = CHEAP_DIST_LIMIT * resScale
        for (i in 0..7) {
            arrayPool.put(i, LinkedList<Array<PointF>>())
        }
    }

    private fun getArray(size: Int) : Array<PointF> {
        return if (USE_POOL) {
            (arrayPool[size] ?: LinkedList<Array<PointF>>().apply {
                arrayPool.put(
                    size,
                    this
                )
            }).poll() ?: Array(size) { PointF() }
        } else {
            Array(size) { PointF() }
        }
    }

    private fun returnArray(arr: Array<PointF>) {
        if (USE_POOL) {
            arrayPool[arr.size].offer(arr)
        }
    }

    fun copy(src: Array<PointF>, dst: Array<PointF>, srcOffset: Int, dstOffset: Int = 0, size: Int) {
        for (i in 0 until size) {
            dst[dstOffset + i] = src[srcOffset + i]
        }
    }

    companion object {
        fun tValue2Scalar(t: Int): Float {
            return t.toFloat() / kMaxTValue
        }

        private const val USE_POOL = false
    }

    data class Segment(
        var distance: Float,
        var ptIndex: Int,
        var tValue: Int = 30,
        var type: SegType = SegType.QUAD
    ) {
        fun getT(): Float {
            return tValue2Scalar(this.tValue)
        }
    }

    data class QuadCoeff(var A: PointF, var B: PointF, var C: PointF) {
        private val fC: PointF
        private val fB: PointF
        private val fA: PointF
        init {
            fC = A
            fB = PointF(2.0f * (B.x - fC.x), 2.0f * (B.y - fC.x))
            fA = PointF(C.x - 2.0f * B.x + fC.x, C.y - 2.0f * B.y + fC.y)
        }
        fun eval(tt: PointF): PointF {
            return PointF(
                (fA.x * tt.x + fB.x) * tt.x + fC.x,
                (fA.y * tt.y + fB.y) * tt.y + fC.y
            )
        }
    }

    private fun quadCoeff(src: Array<PointF>): QuadCoeff {
        return QuadCoeff(src[0], src[1], src[2])
    }

    data class CubicCoeff(var A: PointF, var B: PointF, var C: PointF, var D: PointF) {
        private val fA: PointF
        private val fB: PointF
        private val fC: PointF
        private val fD: PointF
        init {
            fA = PointF(D.x + 3f * (B.x - C.x) - A.x, D.y + 3f * (B.y - C.y) - A.y)
            fB = PointF(3f * (C.x - 2f * B.x + A.x), 3f * (C.y - 2f * B.y + A.y))
            fC = PointF(3f * (B.x - A.x), 3f * (B.y - A.y))
            fD = A
        }
        fun eval(t: Float): PointF {
            return PointF(
                ((fA.x * t + fB.x) * t + fC.x) * t + fD.x,
                ((fA.y * t + fB.y) * t + fC.y) * t + fD.y
            )
        }
    }

    private fun cubicCoeff(src: Array<PointF>): CubicCoeff {
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

    private fun Vec2(point: PointF): Vec2 {
        return Vec2(point.x, point.y)
    }

    enum class SegType {
        LINE,
        QUAD,
        CUBIC,
        CONIC,
    }

    private fun lerp(a: Float, b: Float, t: Float): Float {
        return when {
            t <= 0.0f -> {
                a
            }
            t >= 1.0f -> {
                b
            }
            else -> {
                a + (b - a) * t
            }
        }
    }

    private fun lerp(a: PointF, b: PointF, t: PointF): PointF {

        val x: Float = when {
            t.x <= 0.0f -> {
                a.x
            }
            t.x >= 1.0f -> {
                b.x
            }
            else -> {
                a.x + (b.x - a.x) * t.x
            }
        }

        val y = when {
            t.y <= 0.0f -> {
                a.y
            }
            t.y >= 1.0f -> {
                b.y
            }
            else -> {
                a.y + (b.y - a.y) * t.y
            }
        }

        return PointF(x, y)
    }

    private fun tspanBigEnough(tspan: Int): Boolean {
        return (tspan shr 10) > 0
    }

    private fun quadTooCurvy(pts: Array<PointF>): Boolean {
        // diff = (a/4 + b/2 + c/4) - (a/2 + c/2)
        // diff = -a/4 + b/2 - c/4
        val dx = 0.5f * pts[1].x -
                0.25f * (pts[0].x + pts[2].x)
        val dy = 0.5f * pts[1].y -
                0.25f * (pts[0].y + pts[2].y)

        val dist = max(abs(dx), abs(dy))
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

    private fun cheapDistExceedsLimits(pt: PointF, x: Float, y: Float): Boolean {
        val dist = max(abs(x - pt.x), abs(y - pt.y))
        return dist > tolerance
    }

    private fun validUnitDivide(numerator: Float, denominator: Float): Float {
        var numer = numerator
        var denom = denominator
        if (numer < 0F) {
            numer = -numer
            denom = -denom
        }

        if (denom == 0F || numer == 0F || numer >= denom) {
            return 0F
        }

        val r = numer / denom

        if (r.isNaN()) {
            return 0F
        }

        if (r == 0F) {
            return 0F
        }

        return r
    }

    private fun findQuadMaxCurvature(src: Array<PointF>): Float {
        val ax = src[1].x - src[0].x
        val ay = src[1].y - src[0].y
        val bx = src[0].x - src[1].x - src[1].x + src[2].x
        val by = src[0].y - src[1].y - src[1].y + src[2].y
        return validUnitDivide(-(ax * bx + ay * by), bx * bx + by * by)
    }

    private fun cubicTooCurvy(pts: Array<PointF>): Boolean {
        return cheapDistExceedsLimits(
            pts[1],
            lerp(pts[0].x, pts[3].x, 1.0f / 3.0f),
            lerp(pts[0].y, pts[3].y, 1.0f / 3.0f)
        ) ||
                cheapDistExceedsLimits(
                    pts[2],
                    lerp(pts[0].x, pts[3].x, 2.0f / 3.0f),
                    lerp(pts[0].y, pts[3].y, 2.0f / 3.0f)
                )
    }

    private fun evalQuadAt(src: Array<PointF>, t: Float): PointF {
        val tt = PointF(t, t)
        return quadCoeff(src).eval(tt)
    }

    private fun evalQuadAt(src: Array<PointF>, t: Float, pt: PointF, tangent: Vec2?) {
        pt.set(evalQuadAt(src, t))
        tangent?.set(evalQuadTangentAt(src, t))
    }

    private fun evalQuadTangentAt(src: Array<PointF>, t: Float): Vec2 {
        val p0 = src[0]
        val p1 = src[1]
        val p2 = src[2]

        val b = PointF(
            p1.x - p0.x,
            p1.y - p1.y
        )
        val a = PointF(
            p2.x - p1.x - b.x,
            p2.y - p1.y - b.y
        )

        return Vec2(a.x * t + b.x, a.y * t + b.y)
    }

    private fun quadFoldedLen(pts: Array<PointF>): Float {
        val t = findQuadMaxCurvature(pts)
        val pt = evalQuadAt(pts, t)
        val a = PointF(
            pts[2].x - pt.x,
            pts[2].y - pt.y
        )
        var result = a.length()
        if (0F != t) {
            val b = PointF(
                pts[0].x - pt.x,
                pts[0].y - pt.y
            )
            result += b.length()
        }
        return result
    }

    @Suppress("LocalVariableName")
    fun computeQuadLen(pts: Array<PointF>): Float {
        val a = PointF(0F, 0F)
        val b = PointF(0F, 0F)
        a.x = pts[0].x - 2F * pts[1].x + pts[2].x
        a.y = pts[0].y - 2F * pts[1].y + pts[2].y
        val A = 4F * (a.x * a.x + a.y * a.y)
        if (A == 0F) {
            a.x = pts[2].x - pts[0].x
            a.y = pts[2].y - pts[0].y
            return a.length()
        }
        b.x = 2F * (pts[1].x - pts[0].x)
        b.y = 2F * (pts[1].y - pts[0].y)
        val B = 4F * (a.x * b.x + a.y * b.y)
        val C = (b.x * b.x + b.y * b.y)
        val Sabc = 2F * sqrt(A + B + C)
        val A_2 = sqrt(A)
        val A_32 = 2F * A * A_2
        val C_2 = 2F * sqrt(C)
        val BA = B / A_2
        if (BA + C_2 == 0F) {
            return quadFoldedLen(pts)
        }
        val J = A_32 * Sabc + A_2 * B * (Sabc - C_2)
        val K = 4F * C * A - B * B
        val L = (2F * A_2 + BA + Sabc) / (BA + C_2)
        if (L <= 0F) {
            return quadFoldedLen(pts)
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

        dst[0] = p0
        dst[1] = p01
        dst[2] = lerp(p01, p12, tt)
        dst[3] = p12
        dst[4] = p2
    }

    private fun chopQuadAtHalf(pts: Array<PointF>, dst: Array<PointF>) {
        chopQuadAt(pts, dst, 0.5f)
    }

    private fun computeQuadSegs(
        pts: Array<PointF>, distance: Float,
        mint: Int, maxt: Int, ptIndex: Int
    ): Float {
        var newDist = distance
        if (tspanBigEnough(maxt - mint) && quadTooCurvy(pts)) {
            val tmp = Array(5) { PointF() }
            val halfT = (mint + maxt) shr 1

            chopQuadAtHalf(pts, tmp)
            newDist = distance

            val firstTmp = arrayOf(tmp[0], tmp[1], tmp[2])
            val secondTmp = arrayOf(tmp[2], tmp[3], tmp[4])
            newDist = computeQuadSegs(firstTmp, newDist, mint, halfT, ptIndex)
            newDist = computeQuadSegs(secondTmp, newDist, halfT, maxt, ptIndex)
        } else {
            val d = PointF.length(pts[0].x - pts[2].x, pts[0].y - pts[2].y)
            val prevD = newDist
            newDist += d
            if (newDist > prevD) {
                segments.add(Segment(newDist, ptIndex, maxt, SegType.QUAD))
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

    private fun computeCubicSegs(
        pts: Array<PointF>, distance: Float, mint: Int,
        maxt: Int, ptIndex: Int
    ): Float {
        var newDist = distance
        if (tspanBigEnough(maxt - mint) && cubicTooCurvy(pts)) {
            val tmp = Array(7) { PointF() }
            val halfT = (mint + maxt) shr 1
            chopCubicAtHalf(pts, tmp)

            val firstTmp = arrayOf(tmp[0], tmp[1], tmp[2], tmp[3])
            val secondTmp = arrayOf(tmp[3], tmp[4], tmp[5], tmp[6])
            newDist = computeCubicSegs(firstTmp, distance, mint, halfT, ptIndex)
            newDist = computeCubicSegs(secondTmp, newDist, halfT, maxt, ptIndex)
        } else {
            val d = PointF.length(pts[0].x - pts[3].x, pts[0].y - pts[3].y)
            val prevD = newDist
            newDist += d
            if (newDist > prevD) {
                segments.add(Segment(newDist, ptIndex, maxt, SegType.CUBIC))
            }
        }
        return newDist
    }

    private fun segTo(
        pts: Array<PointF>,
        segType: SegType,
        startT: Float,
        stopT: Float,
        dst: VectorPath
    ) {
        if (startT == stopT) {
            val lastPt = dst.getPoints().last() as PointF
            dst.lineTo(lastPt.x, lastPt.y)
            return
        }

        val tmp0 = Array(7) { PointF() }
        val tmp1 = Array(7) { PointF() }

//        val conic = Conic(pts[0], pts[2], pts[3], pts[1].x)

        when (segType) {
            // Straight segments
            SegType.LINE ->
                if (stopT == 1F) {
                    dst.lineTo(
                        pts[1].x.toDouble(),
                        pts[1].y.toDouble()
                    )
                } else {
                    val stopTT = PointF(stopT, stopT)
                    val pts01 = lerp(pts[0], pts[1], stopTT)
                    dst.lineTo(
                        pts01.x.toDouble(),
                        pts01.y.toDouble()
                    )
                }

            // Quadratic segments
            SegType.QUAD ->
                if (startT == 0F) {
                    if (stopT == 1F) {
                        dst.quadTo(
                            pts[1].x.toDouble(), pts[1].y.toDouble(),
                            pts[2].x.toDouble(), pts[2].y.toDouble()
                        )
                    } else {
                        chopQuadAt(pts, tmp0, stopT)
                        dst.quadTo(
                            tmp0[1].x.toDouble(), tmp0[1].y.toDouble(),
                            tmp0[2].x.toDouble(), tmp0[2].y.toDouble()
                        )
                    }
                } else {
                    chopQuadAt(pts, tmp0, startT)
                    if (stopT == 1F) {
                        dst.quadTo(
                            tmp0[3].x.toDouble(), tmp0[3].y.toDouble(),
                            tmp0[4].x.toDouble(), tmp0[4].y.toDouble()
                        )
                    } else {
                        val secondTmp = arrayOf(tmp0[2], tmp0[3], tmp0[4])
                        chopQuadAt(secondTmp, tmp1, (stopT - startT) / (1F - startT))
                        dst.quadTo(
                            tmp1[1].x.toDouble(), tmp1[1].y.toDouble(),
                            tmp1[2].x.toDouble(), tmp1[2].y.toDouble()
                        )
                    }
                }

            // Conic segments
            SegType.CONIC -> println("Conic SegType not implemented")
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
            SegType.CUBIC ->
                if (startT == 0F) {
                    if (stopT == 1F) {
                        dst.cubicTo(
                            pts[1].x.toDouble(), pts[1].y.toDouble(),
                            pts[2].x.toDouble(), pts[2].y.toDouble(),
                            pts[3].x.toDouble(), pts[3].y.toDouble()
                        )
                    } else {
                        chopCubicAt(pts, tmp0, stopT)
                        dst.cubicTo(
                            tmp0[1].x.toDouble(), tmp0[1].y.toDouble(),
                            tmp0[2].x.toDouble(), tmp0[2].y.toDouble(),
                            tmp0[3].x.toDouble(), tmp0[3].y.toDouble()
                        )
                    }
                } else {
                    chopCubicAt(pts, tmp0, startT)
                    if (stopT == 1F) {
                        dst.cubicTo(
                            tmp0[4].x.toDouble(), tmp0[4].y.toDouble(),
                            tmp0[5].x.toDouble(), tmp0[5].y.toDouble(),
                            tmp0[6].x.toDouble(), tmp0[6].y.toDouble()
                        )
                    } else {
                        val secondTmp = arrayOf(tmp0[3], tmp0[4], tmp0[5], tmp0[6])
                        chopCubicAt(secondTmp, tmp1, (stopT - startT) / (1F - startT))
                        dst.cubicTo(
                            tmp1[1].x.toDouble(), tmp1[1].y.toDouble(),
                            tmp1[2].x.toDouble(), tmp1[2].y.toDouble(),
                            tmp1[3].x.toDouble(), tmp1[3].y.toDouble()
                        )
                    }
                }

        }

    }

    private fun computePosTan(
        pts: Array<PointF>,
        segType: SegType,
        t: Float,
        pos: PointF,
        tangent: Vec2?
    ) {
        when (segType) {
            SegType.LINE -> {
                val tt = PointF(t, t)
                pos.set(lerp(pts[0], pts[1], tt))
                tangent?.normalize(pts[1].x - pts[0].x, pts[1].y - pts[0].y)
            }
            SegType.QUAD -> {
                evalQuadAt(pts, t, pos, tangent)
                tangent?.normalize(tangent.x, tangent.y)
            }
            SegType.CONIC -> println("Conic SegType not implemented")
//            {
//                Conic(pts[0], pts[2], pts[3], pts[1].x).evalAt(t, pos, tangent)
//                if (tangent.isNotEmpty()) {
//                    tangent[0].normalize(tangent[0].x, tangent[0].y)
//                }
//            }
            SegType.CUBIC -> {
                evalCubicAt(pts, t, pos, tangent) // NOTE: Originally, this receives a nullptr,
                // is there anywhere on code where this final
                // value is passed ???
                tangent?.normalize(tangent.x, tangent.y)
            }
        }
    }

    private fun evalCubicAt(
        src: Array<PointF>, t: Float,
        loc: PointF, tangent: Vec2?
    ) {
        loc.set(cubicCoeff(src).eval(t))
        if ((t == 0F && src[0] == src[1]) || (t == 1F && src[2] == src[3])) {
            if (t == 0F) {
                val vec02 = Vec2(
                    src[2].x - src[0].x,
                    src[2].y - src[0].y
                )
                tangent?.set(vec02)
            } else {
                val vec13 = Vec2(
                    src[3].x - src[1].x,
                    src[3].y - src[1].y
                )
                tangent?.set(vec13)
            }
        } else {
            tangent?.set(evalCubicDerivative(src, t))
        }

        // NOTE: This function has a 4th input pointer called curvature, which is not being used
        // anywhere in the code...
    }

    private fun evalCubicDerivative(src: Array<PointF>, t: Float): Vec2 {
        val coeff = QuadCoeff(PointF(), PointF(), PointF())
        val p0 = src[0]
        val p1 = src[1]
        val p2 = src[2]
        val p3 = src[3]

        coeff.A = PointF(
            p3.x + 3F * (p1.x - p2.x) - p0.x,
            p3.y + 3F * (p1.y - p2.y) - p0.y
        )
        coeff.B = PointF(
            2F * (p2.x - 2F * p1.x + p0.x),
            2F * (p2.y - 2F * p1.y + p0.y)
        )
        coeff.C = PointF(
            p1.x - p0.x,
            p1.y - p0.y
        )

        val tt = PointF(t, t)
        return Vec2(coeff.eval(tt))
    }

    private fun buildSegments() {
        val pts = Array(4) { PointF() }
        var ptIndex = this.firstPtIndex
        var distance = 0F
        var isClosed = this.forceClosed

        var prevPt = PointF()
        segments.clear()
        path.visitCmds(
            moveTo = { x, y ->
                ptIndex += 1
                pts[0] = PointF(x.toFloat(), y.toFloat())
                this.pts.add(pts[0])
                prevPt = pts[0]
            },
            lineTo = { x, y ->
                pts[0] = prevPt
                pts[1] = PointF(x.toFloat(), y.toFloat())
                val d = PointF.length(pts[0].x - pts[1].x, pts[0].y - pts[1].y)
                val prevD = distance
                distance += d
                if (distance > prevD) {
                    segments.add(Segment(distance, ptIndex, kMaxTValue, SegType.LINE))
                    this.pts.add(pts[1])
                    ptIndex++
                }
                prevPt = pts[1]
            },
            quadTo = { x0, y0, x1, y1 ->
                pts[0] = prevPt
                pts[1] = PointF(x0.toFloat(), y0.toFloat())
                pts[2] = PointF(x1.toFloat(), y1.toFloat())
                val prevD = distance
                distance = computeQuadSegs(pts, distance, 0, kMaxTValue, ptIndex)
                if (distance > prevD) {
                    this.pts.add(pts[1])
                    this.pts.add(pts[2])
                    ptIndex += 2
                }
                prevPt = pts[2]
            },
            cubicTo = { x0, y0, x1, y1, x2, y2 ->
                pts[0] = prevPt
                pts[1] = PointF(x0.toFloat(), y0.toFloat())
                pts[2] = PointF(x1.toFloat(), y1.toFloat())
                pts[3] = PointF(x2.toFloat(), y2.toFloat())
                val prevD = distance
                distance = computeCubicSegs(pts, distance, 0, kMaxTValue, ptIndex)
                if (distance > prevD) {
                    this.pts.add(pts[1])
                    this.pts.add(pts[2])
                    this.pts.add(pts[3])
                    ptIndex += 3
                }
                prevPt = pts[3]
            },
            close = {
                isClosed = true
            }
        )

        this.pathLength = distance
        this.isClosed = isClosed
        this.firstPtIndex = ptIndex
    }

    fun getLength(): Float {
        if (this.path.isEmpty()) {
            return 0F
        }
        if (this.pathLength < 0F) {
            buildSegments()
        }
        if (pathLength.isNaN()) {
            this.pathLength = 0F
        }
        return this.pathLength
    }

    fun nextContour(): Boolean {
        this.pathLength = -1F
        return getLength() > 0F
    }

    fun setPath(path: Path, forceClosed: Boolean) {
        this.path = path.path
        this.pathLength = -1F   // signal we need to compute it
        this.forceClosed = forceClosed
        this.firstPtIndex = -1
        this.segments.clear()
        this.pts.clear()
    }

    @Suppress("NAME_SHADOWING")
    fun getSegment(startD: Float, stopD: Float, dst: Path, startWithMoveTo: Boolean): Boolean {
        val length = getLength()

        val startD = max(startD, 0F)
        val stopD = min(stopD, length)

        if (startD > stopD) return false
        if (this.segments.count() == 0) return false

        val p = PointF()

        val startPair = distanceToSegment(startD)
        val stopPair = distanceToSegment(stopD)
        var seg = startPair.first
        var startT: Float = startPair.second
        val stopSeg = stopPair.first
        val stopT: Float = stopPair.second
        val allPts = pts.toTypedArray()

        if (startWithMoveTo) {
            val startIndex = seg.first().ptIndex
            val size = pts.size - startIndex
            val tmpPts = getArray(size)
            copy(allPts, tmpPts, startIndex, 0, size)

            // update this.pts
            computePosTan(tmpPts, seg.first().type, startT, p, Vec2(0F, 0F))
            returnArray(tmpPts)
            dst.moveTo(p.x, p.y)
        }

        if (seg.first().ptIndex == stopSeg.first().ptIndex) {
            val startIndex = seg.first().ptIndex
            val size = pts.size - startIndex
            val tmpPts = getArray(size)
            copy(allPts, tmpPts, startIndex, 0, size)

            segTo(tmpPts, seg.first().type, startT, stopT, dst.path)

            returnArray(tmpPts)
            // update this.pts ?
        } else {
            do {
                val startIndex = seg.first().ptIndex
                val size = pts.size - startIndex
                val tmpPts = getArray(size)
                copy(allPts, tmpPts, startIndex, 0, size)

                segTo(tmpPts, seg.first().type, startT, 1F, dst.path)

                returnArray(tmpPts)
                // update this.pts ?
                seg = nextSegment(seg)
                startT = 0F
            } while (seg.first().ptIndex < stopSeg.first().ptIndex)

            val startIndex = seg.first().ptIndex
            val size = pts.size - startIndex
            val tmpPts = getArray(size)
            copy(allPts, tmpPts, startIndex, 0, size)

            segTo(tmpPts, seg.first().type, 0F, stopT, dst.path)

            returnArray(tmpPts)
        }
        return true
    }

    private fun nextSegment(seg: Array<Segment>): Array<Segment> {
        val ptIndex = seg.first().ptIndex
        var index = 0
        do {
            ++index
        } while (seg[index].ptIndex == ptIndex)
        return seg.copyOfRange(index, seg.size)
    }

    private fun findClosestSegment(segments: Array<Segment>, count: Int, distance: Float): Int {
        if (count <= 0) {
            return 0.inv()
        }

        var lo = 0
        var hi = count - 1

        while (lo < hi) {
            val mid = (hi + lo) shr 1
            if (segments[mid].distance < distance) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }

        if (segments[hi].distance < distance) {
            hi += 1
            hi = hi.inv()
        } else if (distance < segments[hi].distance) {
            hi = hi.inv()
        }
        return hi
    }

    private fun distanceToSegment(distance: Float): Pair<Array<Segment>, Float> {
        var index = findClosestSegment(segments.toTypedArray(), segments.count(), distance)
        index = index xor index.shr(31)
        val tmpSeg = segments.toTypedArray().copyOfRange(index, segments.size)

        var startT = 0F
        var startD = 0F
        if (index > 0) {
            startD = segments[index - 1].distance
            if (segments[index - 1].ptIndex == tmpSeg.first().ptIndex) {
                startT = segments[index - 1].getT()
            }
        }

        val t = startT + (tmpSeg.first().getT() - startT) * (distance - startD) /
                (tmpSeg.first().distance - startD)
        return Pair(tmpSeg, t)
    }

    @Suppress("NAME_SHADOWING")
    fun getPosTan(distance: Float, pos: FloatArray, tangent: Vec2?): Boolean {
        if (path.isEmpty()) return false

        val length = getLength()
        val count = segments.count()

        if (count == 0 || length == 0F) return false

        val distance = min(max(0F, distance), length)
        val pairSegT = distanceToSegment(distance)
        val seg = pairSegT.first
        val t = pairSegT.second

        val tmpPts = pts.toTypedArray().copyOfRange(seg.first().ptIndex, pts.size)
        val posPoint = PointF(pos[0], pos[1])
        computePosTan(tmpPts, seg.first().type, t, posPoint, tangent)
        pos[0] = posPoint.x
        pos[1] = posPoint.y
        return true
    }

}