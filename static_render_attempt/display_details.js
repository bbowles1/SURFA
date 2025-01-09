// UORFDetailsPanel.js
const UORFDetailsPanel = ({ region }) => {
  const panelStyle = {
    position: 'absolute',
    right: '1rem',
    top: '1rem',
    padding: '1rem',
    backgroundColor: 'white',
    borderRadius: '0.5rem',
    boxShadow: '0 4px 6px -1px rgb(0 0 0 / 0.1)',
    width: '16rem',
    fontFamily: 'system-ui, -apple-system, sans-serif'
  };

  const labelStyle = {
    fontWeight: 500,
    color: '#374151'
  };

  const containerStyle = {
    display: 'flex',
    flexDirection: 'column',
    gap: '0.5rem',
    fontSize: '0.875rem'
  };

  const exonContainerStyle = {
    paddingLeft: '1rem',
    marginTop: '0.25rem',
    fontSize: '0.75rem'
  };

  const placeholderStyle = {
    fontSize: '0.875rem',
    color: '#6B7280'
  };

  const renderExons = (exons) => {
    return exons.map((exon, i) => (
      React.createElement('div', { key: i },
        `${i + 1}: ${exon.start} - ${exon.end} (${exon.type})`
      )
    ));
  };

  const renderContent = () => {
    if (!region) {
      return React.createElement('div', { style: placeholderStyle },
        'Hover over a region to see details'
      );
    }

    return React.createElement('div', { style: containerStyle }, [
      React.createElement('div', { key: 'id' }, [
        React.createElement('span', { style: labelStyle }, 'ID: '),
        region.id
      ]),
      React.createElement('div', { key: 'type' }, [
        React.createElement('span', { style: labelStyle }, 'Type: '),
        region.type
      ]),
      region.start_codon && React.createElement('div', { key: 'start-codon' }, [
        React.createElement('span', { style: labelStyle }, 'Start Codon: '),
        region.start_codon
      ]),
      React.createElement('div', { key: 'position' }, [
        React.createElement('span', { style: labelStyle }, 'Position: '),
        `${region.start} - ${region.end}`
      ]),
      React.createElement('div', { key: 'length' }, [
        React.createElement('span', { style: labelStyle }, 'Length: '),
        `${region.end - region.start} bp`
      ]),
      region.exons && React.createElement('div', { key: 'exons' }, [
        React.createElement('span', { style: labelStyle }, 'Exons: '),
        region.exons.length,
        React.createElement('div', { style: exonContainerStyle },
          renderExons(region.exons)
        )
      ])
    ]);
  };

  return React.createElement('div', { style: panelStyle },
    renderContent()
  );
};

export default UORFDetailsPanel;