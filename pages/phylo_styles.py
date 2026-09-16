stylesheet = [
    {
        'selector': '.nonterminal',
        'style': {
            'background-opacity': 0,
            "text-halign": "left",
            "text-valign": "top",
        }
    },
    {
        'selector': 'node',
        'style': {
            'font-size': '50px',
            "text-halign": "left",
            "text-valign": "center"
        }
    },

    {
        'selector': '.support',
        'style': {'background-opacity': 0}
    },
    {
        'selector': 'edge',
        'style': {
            "source-endpoint": "inside-to-node",
            "target-endpoint": "inside-to-node",
        }
    },
    {
        'selector': f'edge.path-highlighted',
        'style': {
            'line-color': 'blue',
            'width': 4
        }
    },

    {
        'selector': '.terminal',
        'style': {
            'label': 'data(name)',
            'width': 10,
            'height': 10,
            "text-valign": "center",
            "text-halign": "right",
            'background-color': '#222222'
        }
    },
    {
        'selector': '.colored-terminal',
        'style': {
            'label': 'data(name)',
            'width': 10,
            'height': 10,
            "text-valign": "center",
            "text-halign": "right",
            'background-color': '#222222',

            'text-background-color': 'data(individualColor)',
            'text-background-opacity': 0.5,

            'text-border-color': 'data(individualColor)',
            'text-border-width': 3,
            'text-border-opacity': 1,
            'text-border-style': 'round',

            'text-background-padding': '5px'
        }
    }

]