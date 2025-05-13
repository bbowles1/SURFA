# Goals

My needs are simple:
1. Serve the static home dashboard
2. Pass users calls from the dashboard to the python dataloader.
3. Return the structured JSON from the dataloader.

I want to preserve everything as-is to allow the user to define their build parameters.

# Progress Notes
1. I have to activate the "d3-project" conda env. 
2. My app is named `app.py` so I can serve it via `flask run --debug`.
3. At this time, page just displays a text "index" and does not render anything, when it should be displaying `base.html`. I've followed the "routing" examples to set up the expected structure.
Next steps: Follow templating steps to render the `base.html` page.