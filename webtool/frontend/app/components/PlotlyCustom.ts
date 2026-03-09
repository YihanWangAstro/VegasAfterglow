"use client";

// Custom Plotly bundle with only the trace types we use.
// Reduces JS payload from ~4.8 MB to ~1.5 MB compared to the full plotly.js bundle.
import Plotly from "plotly.js/lib/core";
import Scatter from "plotly.js/lib/scatter";
import ScatterGL from "plotly.js/lib/scattergl";
import Heatmap from "plotly.js/lib/heatmap";
import createPlotlyComponent from "react-plotly.js/factory";

Plotly.register([Scatter, ScatterGL, Heatmap]);

const Plot = createPlotlyComponent(Plotly);
export default Plot;
