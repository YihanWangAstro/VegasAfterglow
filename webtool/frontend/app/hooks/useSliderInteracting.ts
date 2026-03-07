import { useEffect, useRef, useState } from "react";

export function useSliderInteracting() {
  const [sliderInteracting, setSliderInteracting] = useState(false);
  const sliderInteractingRef = useRef(false);

  useEffect(() => {
    const sliderSelector = ".sb-slider input[type='range']";
    const markInteracting = () => {
      if (sliderInteractingRef.current) return;
      sliderInteractingRef.current = true;
      setSliderInteracting(true);
    };
    const clearInteracting = () => {
      if (!sliderInteractingRef.current) return;
      sliderInteractingRef.current = false;
      setSliderInteracting(false);
    };

    const handlePointerDown = (event: PointerEvent) => {
      const target = event.target as Element | null;
      if (!target?.closest(sliderSelector)) return;
      markInteracting();
    };

    window.addEventListener("pointerdown", handlePointerDown, true);
    window.addEventListener("pointerup", clearInteracting, true);
    window.addEventListener("pointercancel", clearInteracting, true);
    window.addEventListener("blur", clearInteracting);

    return () => {
      window.removeEventListener("pointerdown", handlePointerDown, true);
      window.removeEventListener("pointerup", clearInteracting, true);
      window.removeEventListener("pointercancel", clearInteracting, true);
      window.removeEventListener("blur", clearInteracting);
    };
  }, []);

  return sliderInteracting;
}
