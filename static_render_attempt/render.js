export const render = (svg, width, height, uorfs, events = {}) => {
  const margin = {top: 10, right: 10, bottom: 75, left: 10};
  const rectHeight = ((height - margin.top - margin.bottom) / uorfs.regions.length) - 4;
  
  const xScale = d3.scaleLinear()
      .domain([uorfs.start, uorfs.end])
      .range([margin.left, width - margin.right]);
  
  const yScale = d3.scaleBand()
    .domain(uorfs.regions.map(d => d.id))
    .range([height - margin.bottom, margin.top])
    .padding(0.5);  // Add padding to ensure rectangles fit within bounds
  
  const CODON_COLORS = {
    'ATG': '#1f77b4',
    'CTG': '#ff7f0e',
    'TTG': '#2ca02c',
    'GTG': '#d62728',
    'TGG': '#9467bd',
    'TCG': '#8c564b',
    'UTU': '#e377c2',
    'TTT': '#7f7f7f',
    'TTC': '#bcbd22'
  };

  const CODON_STYLES = {
    'ATG': {
        background: '#1f77b4',
        text: '#000000'
    },
    'CTG': {
        background: '#ff7f0e',
        text: '#000000'
    },
    'TTG': {
      background: '#2ca02c',
      text: '#000000'
    }
    'GTG': {
      background: '#d62728',
      text: '#000000'
    },
    'TGG': {
      background: '#9467bd',
      text: '#000000'
    },
    'TCG': {
      background: '#8c564b',
      text: '#000000'
    },
    'UTU': {
      background: '#e377c2',
      text: '#000000'
    },
    'TTT': {
      background: '#7f7f7f',
      text: '#000000'
    },
    'TTC': {
      background: '#bcbd22',
      text: '#000000'
    }
  };

  // Get unique start codons present in the data
  const presentCodons = [...new Set(uorfs.regions
    .map(region => region.start_codon)
    .filter(codon => codon !== undefined && codon !== null))];
  // determine if codons are defined in CODON_COLORS
  const knownCodons = presentCodons.filter(codon => CODON_COLORS[codon]);
  const unknownCodons = presentCodons.filter(codon => !CODON_COLORS[codon]);
        
  // Default color for unknown codons
  const DEFAULT_CODON_COLOR = '#999999';

  // Modified color scale to use codon colors
  const tColorScale = d => {
    const codon = d.start_codon || 'unknown';
    return CODON_COLORS[codon] || DEFAULT_CODON_COLOR;
  };
          

  // Create pattern definitions for UTR regions
  const defs = d3.select(svg).append("defs");

  const legendColor = "#625377";  

  // Create legend pattern
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
  
  // Create patterns for each region
  uorfs.regions.forEach(region => {
    const color = tColorScale(region);
    
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

  const legendGroup = d3.select(svg).append("g")
    .attr("class", "legend");

  // Function to create a legend section
  const createLegendSection = (title, items, yOffset) => {

    // Only create section if there are items to show
    if (items.length === 0) return 0;

    const section = legendGroup.append("g")
      .attr("transform", `translate(0, ${yOffset})`);

    section.append("text")
      .attr("x", 0)
      .attr("y", 0)
      .attr("font-family", "helvetica neue, helvetica, sans-serif")
      .attr("font-size", "12px")
      .attr("font-weight", "bold")
      .text(title);

    const legendSpacing = 20;
    const legendRectWidth = 30;
    const legendRectHeight = 15;

    const legendItem = section.selectAll(".legend-item")
    .data(items)
    .enter()
    .append("g")
    .attr("class", "legend-item")
    .attr("transform", (d, i) => `translate(0, ${i * legendSpacing + 15})`);

  legendItem.append("rect")
    .attr("width", legendRectWidth)
    .attr("height", legendRectHeight)
    .attr("fill", d => d.pattern || d.color)
    .attr("stroke", d => d.stroke || d.color)
    .attr("stroke-width", 1);

  legendItem.append("text")
    .attr("x", legendRectWidth + 5)
    .attr("y", legendRectHeight / 2)
    .attr("dy", "0.35em")
    .attr("font-family", "helvetica neue, helvetica, sans-serif")
    .attr("font-size", "12px")
    .text(d => d.type);

  return section.node().getBBox().height;
  };

// Create region types section
const regionItems = [
  { type: "UTR Region", pattern: "url(#legend-utr-pattern)", stroke: legendColor },
  { type: "CDS Region", pattern: legendColor, stroke: legendColor }
];
const regionHeight = createLegendSection("Region Types", regionItems, 0);

  // Create start codons section for all present codons
  const codonItems = [
    // Known codons with their defined colors
    ...knownCodons.map(codon => ({
      type: codon,
      color: CODON_COLORS[codon]
    })),
    // Unknown codons with default color
    ...unknownCodons.map(codon => ({
      type: codon,
      color: DEFAULT_CODON_COLOR
    }))
  ];

// Only create codon section if there are codons to show
const codonHeight = codonItems.length > 0 ? 
  createLegendSection("Start Codons", codonItems, regionHeight + 20) : 0;

// Get the bounding box of the entire legend
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

// Position the entire legend group in the top-right
legendGroup.attr("transform", `translate(
  ${width - margin.right - legendBBox.width - 10},
  ${margin.top + 10}
)`);





const xAxis = d3.axisBottom(xScale).ticks(5).tickSizeOuter(0);
  
  const chartContainer = d3.select(svg).append('g')
    .attr('class', 'chart-container');
  
  d3.select(svg).append("g")
    .attr("transform", `translate(0,${height - margin.bottom})`)
    .call(xAxis).selectAll("text")
    .style("font-size", "14px");

  // x-axis label
  d3.select(svg).append("text")
    .attr("class", "x-axis-label")
    .attr("text-anchor", "middle")
    .attr("x", width / 2)
    .attr("y", height - margin.bottom / 3)
    .attr("font-family", "helvetica neue, helvetica, sans-serif")
    .attr("font-size", "18px")
    .text("Nucleotide Position");
    
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
        .duration(200)
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
            const color = d3.color(legendColor);
            color.opacity = 0.7;
            return color;
          }
          return tColorScale(regionD);
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
        .delay(200)
        .duration(200)
        .ease(d3.easeLinear)
        .attr('opacity', 1);
    });
};