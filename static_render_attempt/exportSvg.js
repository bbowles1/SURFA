export const initSvgExport = (containerId, filename = 'visualization') => {
  // Add export button to container
  const container = document.getElementById(containerId);
  const exportButton = document.createElement('button');
  exportButton.textContent = 'Export SVG';
  exportButton.className = 'export-button';
  container.appendChild(exportButton);

  const exportSvg = () => {
    // Get the SVG element
    const svgElement = container.querySelector('svg');
    if (!svgElement) {
      console.error('No SVG element found');
      return;
    }

    // Clone the SVG
    const clonedSvg = svgElement.cloneNode(true);
    
    // Ensure SVG has necessary attributes
    clonedSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    clonedSvg.setAttribute('version', '1.1');

    // Define D3's default axis styles
    const d3AxisStyles = `
      .domain {
        fill: none;
        stroke: #000;
      }
      .tick line {
        fill: none;
        stroke: #000;
      }
      .tick text {
        fill: #000;
        font-size: 14px;
        font-family: "helvetica neue", helvetica, sans-serif;
      }
    `;
    
    // Create a style element and add it to the SVG
    const styleElement = document.createElementNS('http://www.w3.org/2000/svg', 'style');
    styleElement.textContent = d3AxisStyles;
    clonedSvg.insertBefore(styleElement, clonedSvg.firstChild);

    // Convert to string with XML declaration
    const serializer = new XMLSerializer();
    const svgString = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' +
                     '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n' +
                     serializer.serializeToString(clonedSvg);

    // Create blob and download link
    const blob = new Blob([svgString], { 
      type: 'image/svg+xml;charset=utf-8'
    });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = `${filename}.svg`;
    
    // Trigger download
    document.body.appendChild(link);
    link.click();
    
    // Cleanup
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  // Add click handler
  exportButton.addEventListener('click', exportSvg);

  // Return cleanup function
  return () => {
    exportButton.removeEventListener('click', exportSvg);
    container.removeChild(exportButton);
  };
};