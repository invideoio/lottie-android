package com.airbnb.lottie.ivvectoranimation

import android.graphics.Path as APath

class PathUtils {
    companion object {
        @JvmStatic
        fun androidPath(path: Path): APath {
            val aPath = APath()
            path.path.visitCmds(
                    moveTo = { x, y ->
                        aPath.moveTo(x.toFloat(), y.toFloat())
                    },
                    lineTo = { x, y ->
                        aPath.lineTo(x.toFloat(), y.toFloat())
                    },
                    quadTo = { x0, y0, x1, y1 ->
                        aPath.quadTo(x0.toFloat(), y0.toFloat(), x1.toFloat(), y1.toFloat())
                    },
                    cubicTo = { x0, y0, x1, y1, x2, y2 ->
                        aPath.cubicTo(x0.toFloat(), y0.toFloat(), x1.toFloat(), y1.toFloat(), x2.toFloat(), y2.toFloat())
                    },
                    close = {
                        aPath.close()
                    }
            )
            return aPath
        }
    }
}