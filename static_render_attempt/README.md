# README

I made this page as a way to preview changes to the uORF rendering function in a way that's version-controlled.

To set up a local server via npm, I used:
1. `npm init` to initialize a project (which created package.json).
2. `npm install vite --save-dev` to install vite.
3. `npm run dev` to deploy a local webserver.


# WIP:

Supposedly this code will apply a cross-hatching pattern to SVG elements:
```
// Create an SVG element
const svg = d3.select("body")
    .append("svg")
    .attr("width", 500)
    .attr("height", 300);

// Define a pattern in the <defs> section
const defs = svg.append("defs");

defs.append("pattern")
    .attr("id", "diagonalHatch")
    .attr("width", 10)
    .attr("height", 10)
    .attr("patternUnits", "userSpaceOnUse")
    .append("path")
    .attr("d", "M0,0 L10,10 M10,0 L0,10")
    .attr("stroke", "#000")
    .attr("stroke-width", 1);

// Draw rectangles with the pattern fill
svg.selectAll("rect")
    .data([50, 120, 190])
    .enter()
    .append("rect")
    .attr("x", (d, i) => i * 120 + 10)
    .attr("y", 50)
    .attr("width", 100)
    .attr("height", 200)
    .attr("fill", "url(#diagonalHatch)");

```