from flask import Flask, send_file
import os

# change dir
os.chdir('/app/static_render_attempt/')

app = Flask(__name__)

@app.route('/')
def serve_index():
    return send_file('/app/static_render_attempt/index.html')

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=True)
