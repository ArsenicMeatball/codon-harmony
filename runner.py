import codon_harmony
from codon_harmony.codon_harmony import AAForm

test = {}
test['input'] = 'test.txt'
form = AAForm(test)
codon_harmony.runner(form)
