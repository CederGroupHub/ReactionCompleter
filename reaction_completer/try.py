import pandas
import pickle

from tqdm import tqdm

from reaction_completer import balance_recipe, CannotBalance

data = pickle.load(open('solgel_data_v9.pkl', 'rb'))
reactions = []

# The Old verison is 96d66414b9b673c46e617e54a1f3e8cd4704c190
tuples = []
for i in tqdm(data):
    for j in i['extracted_data']:
        precursors = j['precursors']
        targets = j['target']

        if not precursors or not targets:
            continue

        text = j['text']
        try:
            reactions = balance_recipe(precursors, targets, text)
            reactions = [x[3] for x in reactions]
        except CannotBalance as e:
            print(e)

        precursor_strings = j['precursors_string']
        target_strings = j['target_string']

        if len(reactions) > 1:
            for r in reactions:
                tuples.append({
                    'precursors': ', '.join(precursor_strings),
                    'targets': ', '.join(target_strings),
                    'reaction': r
                })
        else:
            tuples.append({
                'precursors': ', '.join(precursor_strings),
                'targets': ', '.join(target_strings),
                'reaction': None
            })

data = pandas.DataFrame(tuples, columns=['reaction', 'precursors', 'targets'])
data.to_excel('data.xlsx')
