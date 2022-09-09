# -*- coding: utf-8 -*-
# ************************************************************ #
#   FileName      : web.py
#   Author        : Mengchen
#   Email         : mengchenpu@gmail.com
#   Create on     : 03-11-2022
#   Last modified : 03-11-2022 10:48
#   Version       : V1.0
#   Description   : 
# ************************************************************ #

import os
import time
import subprocess
from flask import Flask, request, render_template

app = Flask(__name__, template_folder="../Code/condensateNetV2/templates")

@app.route("/llpswise/<id>")                       
def unionPriK(id):                               
    return render_template('newnetworkpriC3-0616'+id+'union.html')

if __name__=="__main__":
    app.run(host="0.0.0.0", debug = True, port = 8000)
