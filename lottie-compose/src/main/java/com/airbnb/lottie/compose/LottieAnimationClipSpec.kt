package com.airbnb.lottie.compose

import com.airbnb.lottie.LottieComposition
import com.airbnb.lottie.LottieDrawable

/**
 * Use subclasses of [LottieAnimationClipSpec] to set min/max bounds on the animation playback.
 */
sealed class LottieAnimationClipSpec {

    internal abstract fun applyTo(drawable: LottieDrawable)

    internal abstract fun getMinProgress(composition: LottieComposition): Float

    internal abstract fun getMaxProgress(composition: LottieComposition): Float

    /**
     * Play the animation starting from this frame.
     */
    data class MinFrame(val minFrame: Int) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinFrame(minFrame)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return (minFrame / composition.endFrame).coerceIn(0f, 1f)
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return 1f
        }
    }

    /**
     * Play the animation until this frame.
     */
    data class MaxFrame(val maxFrame: Int, val inclusive: Boolean = true) : LottieAnimationClipSpec() {

        private val actualMaxFrame = if (inclusive) maxFrame else maxFrame - 1

        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMaxFrame(actualMaxFrame)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return 0f
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return (actualMaxFrame / composition.endFrame).coerceIn(0f, 1f)
        }
    }

    data class MinAndMaxFrame(val minFrame: Int, val maxFrame: Int, val maxFrameInclusive: Boolean = true) : LottieAnimationClipSpec() {

        private val actualMaxFrame = if (maxFrameInclusive) maxFrame else maxFrame - 1

        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinAndMaxFrame(minFrame, actualMaxFrame)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return (minFrame / composition.endFrame).coerceIn(0f, 1f)
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return (actualMaxFrame / composition.endFrame).coerceIn(0f, 1f)
        }
    }

    data class MinProgress(val minProgress: Float) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinProgress(minProgress)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return minProgress
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return 1f
        }
    }

    data class MaxProgress(val maxProgress: Float) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMaxProgress(maxProgress)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return 0f
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return maxProgress
        }
    }

    data class MinAndMaxProgress(val minProgress: Float, val maxProgress: Float) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinAndMaxProgress(minProgress, maxProgress)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return minProgress
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return maxProgress
        }
    }

    data class MinMarker(val minMarker: String) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinFrame(minMarker)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return ((composition.getMarker(minMarker)?.startFrame ?: 0f) / composition.endFrame).coerceIn(0f, 1f)
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return 1f
        }
    }

    data class MaxMarker(val maxMarker: String) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMaxFrame(maxMarker)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return 0f
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return ((composition.getMarker(maxMarker)?.startFrame ?: 0f) / composition.endFrame).coerceIn(0f, 1f)
        }
    }

    data class MinAndMaxMarker(val minMarker: String, val maxMarker: String, val playMaxMarkerStartFrame: Boolean = true) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinAndMaxFrame(minMarker, maxMarker, playMaxMarkerStartFrame)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return ((composition.getMarker(minMarker)?.startFrame ?: 0f) / composition.endFrame).coerceIn(0f, 1f)
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            return ((composition.getMarker(maxMarker)?.startFrame ?: 0f) / composition.endFrame).coerceIn(0f, 1f)
        }
    }

    data class Marker(val marker: String) : LottieAnimationClipSpec() {
        override fun applyTo(drawable: LottieDrawable) {
            drawable.setMinAndMaxFrame(marker)
        }

        override fun getMinProgress(composition: LottieComposition): Float {
            return ((composition.getMarker(marker)?.startFrame ?: 0f) / composition.endFrame).coerceIn(0f, 1f)
        }

        override fun getMaxProgress(composition: LottieComposition): Float {
            val marker = composition.getMarker(marker) ?: return 1f
            return ((marker.startFrame + marker.durationFrames) / composition.endFrame).coerceIn(0f, 1f)
        }
    }
}