from flask import Flask
from flask import url_for

app = Flask(__name__)

@app.route('/')
def index():
    return 'index'

@app.route('/uorfs/gene_id')
def uorf_view(gene_id):
    return f'uORF view for {gene_id}.'

with app.test_request_context():
    print(url_for('index'))
    print(url_for('uorf_view', gene_id='ENSG00000081189'))