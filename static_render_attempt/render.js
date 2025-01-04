export const render = (svg, width, height, uorfs) => {
  const margin = ({top: 10, right: 10, bottom: 25, left: 10});
  const rectHeight = ((height - margin.top - margin.bottom) / uorfs.regions.length ) - 4;
  
  const xScale = d3.scaleLinear()
      .domain([uorfs.start, uorfs.end])
      .range([margin.left, width - margin.right]);
  
  const yScale = d3.scaleBand()
      .domain(uorfs.regions.map(d => d.id))
      .range([height - margin.bottom, margin.top]);
  
  const tColorScale = d3.scaleOrdinal()
      .domain(uorfs.regions.map(d => d.id))
      .range(['#855C75', '#D9AF6B','#AF6458','#736F4C','#526A83','#625377','#68855C','#9C9C5E','#A06177','#8C785D']);
  
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
        // Forward diagonal lines
        g.append("path").attr("d", "M-2,2 l4,-4 M0,8 l8,-8 M6,10 l4,-4");
        // Backward diagonal lines
        g.append("path").attr("d", "M-2,6 l8,-8 M0,8 l8,-8 M6,2 l4,-4");
      });
  });
  
  const xAxis = d3.axisBottom(xScale).ticks(5).tickSizeOuter(0);
  
  const chartContainer = d3.select(svg).append('g').attr('class', 'chart-container');
  
  d3.select(svg).append("g")
      .attr("transform", `translate(0,${600 - 20})`)
      .call(xAxis);
  
  const tooltip = d3.select(svg).append('text').attr('class', 'tooltip-container')
      .attr('dy', -10)
      .style("font-family", "helvetica neue, helvetica, sans-serif")
      .style("font-size", ".8rem")
      .style("fill", "#666")
      .style('display', 'none');
  
  chartContainer.selectAll('g')
    .data(uorfs.regions)
    .enter().append('g').append('line')
        .attr('x1', d => xScale(d.start))
        .attr('x2', d => xScale(d.start))
        .attr('y1', d => yScale(d.id))
        .attr('y2', d => yScale(d.id))
        .attr('stroke-width', 1)
        .attr('stroke', d => tColorScale(d.id))
        .transition().duration(400)
            .ease(d3.easeQuad)
            .attr('x2', d => xScale(d.end))
        .each(function (regionD) {
            const node = d3.select(this.parentNode);
            node.selectAll('rect')
              .data(d => d.exons)
              .enter().append('rect')
                  .attr('opacity', 0)
                  .attr('stroke-width', 1)
                  .attr('width', d => xScale(d.end) - xScale(d.start))
                  .attr('height', rectHeight)
                  .attr('x', d => xScale(d.start)) 
                  .attr('y', yScale(regionD.id) - rectHeight / 2)
                  .attr('fill', d => {
                    if (d.type === 'utr') {
                      // Use cross-hatch pattern for UTR
                      return `url(#cross-hatch-${regionD.id})`;
                    } else if (d.type === 'cds') {
                      // Semi-transparent solid color for CDS
                      const color = d3.color(tColorScale(regionD.id));
                      return color.copy({opacity: 0.7});
                    }
                    // Solid color for uorf
                    return tColorScale(regionD.id);
                  })
                  .attr('stroke', d => tColorScale(regionD.id))
                  .on('mouseover', () => { tooltip.style('display', null); })
                  .on('mouseout', () => { tooltip.style('display', 'none'); })
                  .on('mousemove', () => {
                    const mouseX = d3.pointer(event)[0];
                    const mouseY = d3.pointer(event)[1];
                    tooltip.attr('transform', `translate(${mouseX + 5}, ${mouseY + 5})`);
                    tooltip.text(`Transcript id: ${regionD.id}`)
                    if (mouseX > width / 2) {
                        tooltip.attr('transform', `translate(${mouseX + 10}, ${mouseY + 5})`)
                          .attr('text-anchor', 'end');
                    } else {
                        tooltip.attr('transform', `translate(${mouseX + 15}, ${mouseY + 5})`)
                          .attr('text-anchor', 'start');
                    }
                    if (mouseY > height / 2) {
                        tooltip.attr('transform', `translate(${mouseX + 5}, ${mouseY + 5})`)
                    } else {
                        tooltip.attr('transform', `translate(${mouseX - 10}, ${mouseY + 25})`)
                    }
                  })
                  .transition().delay(1500).duration(1200)
                      .ease(d3.easeLinear)
                      .attr('opacity', 1)
        });
}