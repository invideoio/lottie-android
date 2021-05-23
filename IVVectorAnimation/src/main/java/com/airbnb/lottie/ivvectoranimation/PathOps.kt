package com.airbnb.lottie.ivvectoranimation

class PathOps {
/*fun OpDebug(one: Path, two: Path, ops: Array<SkPathOp>, result: Path) : Boolean {
    SkSTArenaAlloc<4096> allocator;  // FIXME: add a constant expression here, tune
    SkOpContour contour;
    SkOpContourHead* contourList = static_cast<SkOpContourHead*>(&contour);
    SkOpGlobalState globalState(contourList, &allocator
            SkDEBUGPARAMS(skipAssert) SkDEBUGPARAMS(testName));
    SkOpCoincidence coincidence(&globalState);

    var op = ops[0]
    val opIndex = op.ordinal

    op = gOpInverse[opIndex][if (one.isInverseFillType()) 1 else 0][if (two.isInverseFillType()) 1 else 0]
    val inverseFill = gOutInverse[opIndex][if (one.isInverseFillType()) 1 else 0][if (two.isInverseFillType()) 1 else 0]
    val fillType = if (inverseFill)  Path.FillType.INVERSE_EVEN_ODD else Path.FillType.EVEN_ODD

//    SkScalar scaleFactor = SkTMax(ScaleFactor(one), ScaleFactor(two));
//    SkPath scaledOne, scaledTwo;
    var minuend = one
    var subtrahend = two
*//*    if (scaleFactor > SK_Scalar1) {
        ScalePath(one, 1.f / scaleFactor, &scaledOne);
        minuend = &scaledOne;
        ScalePath(two, 1.f / scaleFactor, &scaledTwo);
        subtrahend = &scaledTwo;
    } else {
        minuend = &one;
        subtrahend = &two;
    }*//*
    if (op == SkPathOp.kReverseDifference_SkPathOp) {
        val tmp = minuend
        minuend = subtrahend
        subtrahend = tmp
        op = SkPathOp.kDifference_SkPathOp;
    }

    // turn path into list of segments
    SkOpEdgeBuilder builder(*minuend, contourList, &globalState);
    if (builder.unparseable()) {
        return false;
    }
    const int xorMask = builder.xorMask();
    builder.addOperand(*subtrahend);
    if (!builder.finish()) {
        return false;
    }

    const int xorOpMask = builder.xorMask();
    if (!SortContourList(&contourList, xorMask == kEvenOdd_PathOpsMask,
            xorOpMask == kEvenOdd_PathOpsMask)) {
        result->reset();
        result->setFillType(fillType);
        return true;
    }
    // find all intersections between segments
    SkOpContour* current = contourList;
    do {
        SkOpContour* next = current;
        while (AddIntersectTs(current, next, &coincidence)
                && (next = next->next()))
            ;
    } while ((current = current->next()));

    bool success = HandleCoincidence(contourList, &coincidence);

    if (!success) {
        return false;
    }

    // construct closed contours
    result->reset();
    result->setFillType(fillType);
    SkPathWriter wrapper(*result);
    if (!bridgeOp(contourList, op, xorMask, xorOpMask, &wrapper)) {
        return false;
    }
    wrapper.assemble();  // if some edges could not be resolved, assemble remaining

    if (scaleFactor > 1) {
        ScalePath(*result, scaleFactor, result);
    }
    return true;
}*/

/*    fun ScaleFactor(path: Path) {
        var largest = 0

        val oneBounds = path.bounds.left
        for (index in 0..4) {
            largest = max(largest, abs(oneBounds[index]));
        }
    for (int index = 0; index < 4; ++index) {
        largest = SkTMax(largest, SkScalarAbs(oneBounds[index]));
    }
    SkScalar scale = twoTo10;
    SkScalar next;
    while ((next = scale * twoTo10) < largest) {
        scale = next;
    }
    return scale == twoTo10 ? SK_Scalar1 : scale;
}*/

    // diagram of why this simplifcation is possible is here:
// https://skia.org/dev/present/pathops link at bottom of the page
// https://drive.google.com/file/d/0BwoLUwz9PYkHLWpsaXd0UDdaN00/view?usp=sharing
    private val gOpInverse =
            arrayOf(
                    arrayOf(
                            arrayOf(SkPathOp.kDifference_SkPathOp, SkPathOp.kIntersect_SkPathOp),
                            arrayOf(SkPathOp.kUnion_SkPathOp, SkPathOp.kReverseDifference_SkPathOp)
                    ),
                    arrayOf(
                            arrayOf(SkPathOp.kIntersect_SkPathOp, SkPathOp.kDifference_SkPathOp),
                            arrayOf(SkPathOp.kReverseDifference_SkPathOp, SkPathOp.kUnion_SkPathOp)
                    ),
                    arrayOf(
                            arrayOf(SkPathOp.kUnion_SkPathOp, SkPathOp.kReverseDifference_SkPathOp),
                            arrayOf(SkPathOp.kDifference_SkPathOp, SkPathOp.kIntersect_SkPathOp)
                    ),
                    arrayOf(
                            arrayOf(SkPathOp.kXOR_SkPathOp, SkPathOp.kXOR_SkPathOp),
                            arrayOf(SkPathOp.kXOR_SkPathOp, SkPathOp.kXOR_SkPathOp)
                    ),
                    arrayOf(
                            arrayOf(SkPathOp.kReverseDifference_SkPathOp, SkPathOp.kUnion_SkPathOp),
                            arrayOf(SkPathOp.kIntersect_SkPathOp, SkPathOp.kDifference_SkPathOp)
                    )
            )

    private val gOutInverse =
            arrayOf(
                    arrayOf(
                            arrayOf(false, false),
                            arrayOf(true, false)
                    ),
                    arrayOf(
                            arrayOf(false, false),
                            arrayOf(false, true)
                    ),
                    arrayOf(
                            arrayOf(false, true),
                            arrayOf(true, true)
                    ),
                    arrayOf(
                            arrayOf(false, true),
                            arrayOf(true, false)
                    ),
                    arrayOf(
                            arrayOf(false, true),
                            arrayOf(false, false)
                    )
            )

    enum class SkPathOp {
        kDifference_SkPathOp,         //!< subtract the op path from the first path
        kIntersect_SkPathOp,          //!< intersect the two paths
        kUnion_SkPathOp,              //!< union (inclusive-or) the two paths
        kXOR_SkPathOp,                //!< exclusive-or the two paths
        kReverseDifference_SkPathOp,  //!< subtract the first path from the op path
    }

    private fun Path.isInverseFillType(): Boolean {
        return fillType == Path.FillType.INVERSE_EVEN_ODD || fillType == Path.FillType.INVERSE_WINDING
    }

    companion object {
        const val twoTo10 = 1024f
    }
}
