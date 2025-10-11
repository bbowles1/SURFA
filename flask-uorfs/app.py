from flask import Flask, url_for, render_template
from flask import jsonify, send_from_directory # extra things to try and use for POST
import json

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('base.html')

@app.route('/uorfs/<gene_id>')
def uorf_view(gene_id):
    return f"uORF view for {gene_id}."

with app.test_request_context():
    print(url_for('index'))
    print(url_for('uorf_view', gene_id='ENSG00000081189'))

# load JSON
@app.route('/api/uorf')
def uorf():
    with open('data/uorfs.json', 'r') as f:
        return jsonify(json.load(f))


@app.route('/static/<path:path>')
def send_static(path):
    return send_from_directory('static/', path)