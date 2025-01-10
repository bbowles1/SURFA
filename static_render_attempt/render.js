export const render = (svg, width, height, uorfs, events = {}) => {
  const margin = {top: 10, right: 10, bottom: 25, left: 10};
  const rectHeight = ((height - margin.top - margin.bottom) / uorfs.regions.length) - 4;
  
  const xScale = d3.scaleLinear()
      .domain([uorfs.start, uorfs.end])
      .range([margin.left, width - margin.right]);
  
  const yScale = d3.scaleBand()
      .domain(uorfs.regions.map(d => d.id))
      .range([height - margin.bottom, margin.top]);
  
  const tColorScale = d3.scaleOrdinal()
      .domain(uorfs.regions.map(d => d.id))
      .range(['#855C75', '#D9AF6B', '#AF6458', '#736F4C', '#526A83', '#625377', '#68855C', '#9C9C5E', '#A06177', '#8C785D']);


  // Create pattern definitions for UTR regions
  const defs = d3.select(svg).append("defs");
  
  // Create a diagonal line pattern for each possible color
  uorfs.regions.forEach(region => {
    const color = tColorScale(region.id);
    
    // Pattern for forward diagonal lines
    defs.append("pattern")
      .attr("id", `diagonal-${region.id}`)
      .attr("patternUnits", "userSpaceOnUse")
      .attr("width", 8)
      .attr("height", 8)
      .append("path")
      .attr("d", "M-2,2 l4,-4 M0,8 l8,-8 M6,10 l4,-4")
      .attr("stroke", color)
      .attr("stroke-width", 1)
      .attr("stroke-opacity", 0.7);
      
    // Pattern for backward diagonal lines (creates cross-hatch)
    defs.append("pattern")
      .attr("id", `cross-hatch-${region.id}`)
      .attr("patternUnits", "userSpaceOnUse")
      .attr("width", 8)
      .attr("height", 8)
      .append("g")
      .attr("fill", "none")
      .attr("stroke", color)
      .attr("stroke-width", 1)
      .attr("stroke-opacity", 0.7)
      .call(g => {
        g.append("path").attr("d", "M-2,2 l4,-4 M0,8 l8,-8 M6,10 l4,-4");
        g.append("path").attr("d", "M-2,6 l8,-8 M0,8 l8,-8 M6,2 l4,-4");
      });
  });

  
  // Create legend group without translation first
  const legendGroup = d3.select(svg).append("g")
    .attr("class", "legend");

  const legendColor = "#526A83";

  // Add pattern for legend UTR
  defs.append("pattern")
    .attr("id", "legend-utr-pattern")
    .attr("patternUnits", "userSpaceOnUse")
    .attr("width", 8)
    .attr("height", 8)
    .append("g")
    .attr("fill", "none")
    .attr("stroke", legendColor)
    .attr("stroke-width", 1)
    .attr("stroke-opacity", 0.7)
    .call(g => {
      g.append("path").attr("d", "M-2,2 l4,-4 M0,8 l8,-8 M6,10 l4,-4");
      g.append("path").attr("d", "M-2,6 l8,-8 M0,8 l8,-8 M6,2 l4,-4");
    });

  // Add legend title
  legendGroup.append("text")
    .attr("x", 0)
    .attr("y", 0)
    .attr("font-family", "helvetica neue, helvetica, sans-serif")
    .attr("font-size", "12px")
    .attr("font-weight", "bold")
    .text("Region Types");

  const legendItems = [
    { type: "UTR Region", pattern: "url(#legend-utr-pattern)" },
    { type: "CDS Region", pattern: legendColor }
  ];

  const legendSpacing = 20;
  const legendRectWidth = 30;
  const legendRectHeight = 15;

  // Add legend items
  const legendItem = legendGroup.selectAll(".legend-item")
    .data(legendItems)
    .enter()
    .append("g")
    .attr("class", "legend-item")
    .attr("transform", (d, i) => `translate(0, ${i * legendSpacing + 15})`);

  // Add rectangles for each legend item
  legendItem.append("rect")
    .attr("width", legendRectWidth)
    .attr("height", legendRectHeight)
    .attr("fill", d => d.pattern)
    .attr("stroke", legendColor)
    .attr("stroke-width", 1);

  // Add text labels for each legend item
  legendItem.append("text")
    .attr("x", legendRectWidth + 5)
    .attr("y", legendRectHeight / 2)
    .attr("dy", "0.35em")
    .attr("font-family", "helvetica neue, helvetica, sans-serif")
    .attr("font-size", "12px")
    .text(d => d.type);

  // Get the bounding box of the legend
  const legendBBox = legendGroup.node().getBBox();
  
  // Add semi-transparent white background for legend
  legendGroup.insert("rect", ":first-child")
    .attr("x", -5)
    .attr("y", -5)
    .attr("width", legendBBox.width + 10)
    .attr("height", legendBBox.height + 10)
    .attr("fill", "white")
    .attr("fill-opacity", 0.9)
    .attr("rx", 5)
    .attr("ry", 5);

  // Now position the entire legend group in the top-right
  legendGroup.attr("transform", `translate(
    ${width - margin.right - legendBBox.width - 10},
    ${margin.top + 10}
  )`);
  

  const xAxis = d3.axisBottom(xScale).ticks(5).tickSizeOuter(0);
  
  const chartContainer = d3.select(svg).append('g')
    .attr('class', 'chart-container');
  
  d3.select(svg).append("g")
    .attr("transform", `translate(0,${height - margin.bottom})`)
    .call(xAxis);
  
  const tooltip = d3.select(svg).append('text')
    .attr('class', 'tooltip-container')
    .attr('dy', -10)
    .style("font-family", "helvetica neue, helvetica, sans-serif")
    .style("font-size", ".8rem")
    .style("fill", "#666")
    .style('display', 'none');
  
  chartContainer.selectAll('g')
    .data(uorfs.regions)
    .enter()
    .append('g')
    .each(function(regionD) {
      const group = d3.select(this);
      
      // Add the baseline
      group.append('line')
        .attr('x1', d => xScale(d.start))
        .attr('x2', d => xScale(d.start))
        .attr('y1', d => yScale(d.id))
        .attr('y2', d => yScale(d.id))
        .attr('stroke-width', 1)
        .attr('stroke', d => tColorScale(d.id))
        .transition()
        .duration(400)
        .ease(d3.easeQuad)
        .attr('x2', d => xScale(d.end));
      
      // Add the exon rectangles
      group.selectAll('rect')
        .data(d => d.exons)
        .enter()
        .append('rect')
        .attr('opacity', 0)
        .attr('stroke-width', 1)
        .attr('width', d => xScale(d.end) - xScale(d.start))
        .attr('height', rectHeight)
        .attr('x', d => xScale(d.start))
        .attr('y', yScale(regionD.id) - rectHeight / 2)
        .attr('fill', d => {
          if (d.type === 'utr') {
            return `url(#cross-hatch-${regionD.id})`;
          } else if (d.type === 'cds') {
            const color = d3.color(tColorScale(regionD.id));
            color.opacity = 0.7;
            return color;
          }
          return tColorScale(regionD.id);
        })
        .attr('stroke', d => tColorScale(regionD.id))
        .on('mouseover', (event) => {
          tooltip.style('display', null);
          if (events.onRegionHover) {
            events.onRegionHover(regionD);
          }
        })
        .on('mouseout', (event) => {
          tooltip.style('display', 'none');
          if (events.onRegionLeave) {
            events.onRegionLeave();
          }
        })
        .on('mousemove', (event) => {
          const [mouseX, mouseY] = d3.pointer(event);
          let tooltipX = mouseX + 5;
          let tooltipY = mouseY + 5;
          
          if (mouseX > width / 2) {
            tooltipX = mouseX - 10;
            tooltip.attr('text-anchor', 'end');
          } else {
            tooltip.attr('text-anchor', 'start');
          }
          
          if (mouseY > height / 2) {
            tooltipY = mouseY - 10;
          }
          
          tooltip
            .attr('transform', `translate(${tooltipX},${tooltipY})`)
            .text(`Transcript id: ${regionD.id}`);
        })
        .transition()
        .delay(1500)
        .duration(1200)
        .ease(d3.easeLinear)
        .attr('opacity', 1);
    });
};