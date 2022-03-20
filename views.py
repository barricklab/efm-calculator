from django.shortcuts import render
from django.http import HttpResponseRedirect
from django import forms
from django.core.files.base import ContentFile

import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tempfile import NamedTemporaryFile
import os
#from efm_helper import *

from efm.efm_helper import *

# Create your views here.
class SeqForm(forms.Form):
    '''
    This defines the DNA sequence submission form
    '''

    sequence_file = forms.FileField(required=False)
    sequence = forms.CharField(widget=forms.Textarea, required=False, min_length=50, max_length=50000)
    organism = forms.ChoiceField(choices=(
        ('ecoli', 'E. coli'),
        ('reca', 'E. coli, RecA-'),
        # ('yeast', 'S. cerevisae'),
    ))
    check_features = forms.BooleanField(label='Annotated Regions Only', required=False)
    check_features.help_text = 'Only assign rates for hypermutable sites that overlap ' \
                               'with annotated regions. Sequence must be annotated.'
 
    def clean(self):
        '''
        This is where we detect what type of sequence has been uploaded or placed in in the textbox.
        All form validation happens here.
        '''
        cleaned_data = super(SeqForm, self).clean()
        sequence_file = cleaned_data.get("sequence_file")
        sequence = cleaned_data.get("sequence")

        if sequence_file and sequence:
            raise forms.ValidationError('Please either upload a sequence or copy and paste it into the textbox. '
                                        'Do not submit both.')
        elif not (sequence_file or sequence):
            raise forms.ValidationError('Please submit a sequence for analysis.')
        elif sequence_file:
            if sequence_file.size < 50000:
                with NamedTemporaryFile(delete=False) as data_file:
                    for chunk in sequence_file.chunks():
                        data_file.write(chunk)
            else:
                raise forms.ValidationError('File size is too large.')
        elif sequence:
            with NamedTemporaryFile(delete=False) as data_file:
                data_file.write(sequence.strip().encode(encoding='UTF-8'))

        del cleaned_data['sequence']
        del cleaned_data['sequence_file']

        with open(data_file.name, 'rU') as temp_file:
            cleaned_data['temp_file'] = temp_file.name
            lines = temp_file.readlines()
            data = lines[0].strip()

            if re.match(r'^LOCUS', data):
                # Format is genbank...
                try:
                    genome = SeqIO.read(temp_file.name, 'genbank')
                    cleaned_data['features'] = get_genbank_features(genome)
                    cleaned_data['raw_sequence'] = str(genome.seq).upper()
                    cleaned_data['title'] = genome.id
                    output_handle = NamedTemporaryFile(delete=False)
                    SeqIO.write(SeqRecord(Seq(cleaned_data['raw_sequence']),
                                          id=cleaned_data['title']), output_handle, "fasta")
                    output_handle.close()
                    cleaned_data['fasta_file'] = output_handle.name
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: Genbank filetype detected, but file is malformed.')
            elif re.match(r'^\>', data):
                # Format is FASTA...
                try:
                    genome = SeqIO.read(temp_file.name, 'fasta')
                    cleaned_data['raw_sequence'] = str(genome.seq).upper()
                    cleaned_data['title'] = genome.id
                    cleaned_data['fasta_file'] = temp_file.name
                    cleaned_data['features'] = ''
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: FASTA filetype detected, but file is malformed.')

            elif re.match(r"^<\?[x][m][l].+\?>", data):
                # Format is BioBrick XML...
                try:
                    tree = ElementTree.parse(temp_file.name)
                    cleaned_data['features'] = get_biobrick_features(tree)
                    for node in tree.iter(tag='sequences'):
                        cleaned_data['raw_sequence'] = ''.join(node.find('seq_data').text.split()).upper()
                    for node in tree.iter(tag='part'):
                        cleaned_data['title'] = node.find('part_name').text
                    output_handle = NamedTemporaryFile(delete=False)
                    SeqIO.write(SeqRecord(Seq(cleaned_data['raw_sequence']),
                                          id=cleaned_data['title']), output_handle, "fasta")
                    output_handle.close()
                    cleaned_data['fasta_file'] = output_handle.name
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: BioBrick XML file could not be processed.')

            elif re.match(r'^[ATCGatcg]+', data):
                # Format is just raw sequence
                try:
                    joined_lines = ''.join(lines)
                    cleaned_data['raw_sequence'] = ''.join(joined_lines.split()).upper()
                    cleaned_data['title'] = 'Untitled'
                    cleaned_data['features'] = ''
                    output_handle = NamedTemporaryFile(delete=False)
                    SeqIO.write(SeqRecord(Seq(cleaned_data['raw_sequence']),
                                          id=cleaned_data['title']), output_handle, "fasta")
                    output_handle.close()
                    cleaned_data['fasta_file'] = output_handle.name
                except:
                    if os.path.isfile(cleaned_data['temp_file']):
                        os.remove(cleaned_data['temp_file'])
                    raise forms.ValidationError('Error: FASTA file could not be written from textbox input.')

            else:
                if os.path.isfile(cleaned_data['temp_file']):
                    os.remove(cleaned_data['temp_file'])
                raise forms.ValidationError('The submitted sequence is not valid and cannot be processed.')

            if cleaned_data['check_features'] is True and cleaned_data['features'] == '':
                raise forms.ValidationError('The submitted sequence contains no annotations. Please uncheck ' \
                                            '"Annotated Regions Only".')

        return cleaned_data


def get_sequence(request):
    '''
    This view processes the SeqForm defined above.
    '''
    if request.method == 'POST':
        form = SeqForm(request.POST, request.FILES)
        if form.is_valid():
            results = process_efm(form)
            os.remove(form.cleaned_data.get('fasta_file'))
            if os.path.isfile(form.cleaned_data.get('temp_file')):
                os.remove(form.cleaned_data.get('temp_file'))
            return render(request, 'efm/results.html', results)
    else:
        form = SeqForm()

    return render(request, 'efm/form.html', {'form': form, 'version': EFM_VERSION})
