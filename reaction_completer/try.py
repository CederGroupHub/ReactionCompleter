import pickle

from reaction_completer import balance_recipe, CannotBalance

data = pickle.load(open('solgel_data_v9_100.pkl', 'rb'))
reactions = []

for i in data:
    for j in i['extracted_data']:
        precursors = j['precursors']
        targets = j['target']
        text = j['text']
        try:
            reactions = balance_recipe(precursors, targets, text)
        except CannotBalance as e:
            print(e)

        print(reactions, precursors, targets)
