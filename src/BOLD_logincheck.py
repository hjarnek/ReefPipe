import argparse, pkg_resources
from boldigger_cline.login import login

# Create parser
parser = argparse.ArgumentParser(description='Check connection to BOLDsystems.org')

# Add arguments
parser.add_argument('-u', '--username', required = True, help = 'Specify the BOLDsystems.org username.')
parser.add_argument('-p', '--password', required = True, help = 'Specify the BOLDsystems.org password.')

# Parse the arguments
args = parser.parse_args()

def checklogin(username, password):
    certs = pkg_resources.resource_filename('boldigger_cline.login', 'data/certs.pem')
    session = login(username, password, certs)
    return session

checklogin(args.username, args.password)
