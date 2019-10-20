"""
submit requests made in LCOGT_make_19B20A_requests.py
"""
##########
import pickle, requests, socket

if socket.gethostname() == 'brik':
    api_file = '/home/luke/.lcogt_api_token'
elif 'astro' in socket.gethostname():
    api_file = '/Users/luke/.lcogt_api_token'
else:
    raise NotImplementedError('where to get API file?')

with open(api_file, 'r') as f:
    l = f.readlines()
token = str(l[0].replace('\n',''))
##########

pkl_savpath = (
    '../results/LCOGT_19B20A_observability/all_requests_19B.pkl'
)
with open(pkl_savpath, 'rb') as f:
    r = pickle.load(f)

# DEMO WOULD LOOK LIKE THIS FIXME
requestgroup = r[0][0]

# Submit the fully formed RequestGroup
response = requests.post(
    'https://observe.lco.global/api/requestgroups/',
    headers={'Authorization': 'Token {}'.format(token)},
    json=requestgroup  # Make sure you use json!
)

# Make sure the API call was successful
try:
    response.raise_for_status()
except requests.exceptions.HTTPError as exc:
    print('API call failed: {}'.format(response.content))
    raise exc

requestgroup_dict = response.json()  # The API will return the newly submitted requestgroup as json

# Print out the url on the portal where we can view the submitted request
print('View the observing request: https://observe.lco.global/requestgroups/{}/'.format(requestgroup_dict['id']))
