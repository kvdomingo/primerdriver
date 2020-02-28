/* jshint esversion: 6 */

document.addEventListener('DOMContentLoaded', () => {

    $('#sequence').on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        var cursorposition = $('#sequence').prop('selectionStart');
        document.getElementById('cursorposition').innerHTML = cursorposition + 1;
        document.getElementById('counter').innerHTML = `${this.value.length}/8000`;
    });

    $('#form').on('keyup keypress keydown mousedown mouseup change blur', function(e) {
        if ($('#mut_type').val() === 'sub') {
            $('#target').prop('required', true);
        }
        else {
            $('#target').prop('required', false);
        }

        var valid_sequence = ($('#sequence').val().length > 0);
        var valid_target = ($('#target').val().length > 0);
        var valid_position = !(isNaN(parseInt($('#position').val()))) && (parseInt($('#position').val()) > 0);
        var valid_replacement = ($('#target').val().length > 0);
        var valid_type = ($('#mut_type').val().length > 0);
        if (valid_sequence && valid_target && valid_position && valid_replacement && valid_type) {
            $('#submit').removeClass('disabled');
        }
        else {
            $('#submit').addClass('disabled');
        }
    });

    document.querySelector('#submit').onclick = () => {
        const xhttp = new XMLHttpRequest();
        const sequence = document.querySelector('#sequence').value;
        const target = document.querySelector('#target').value;
        const position = document.querySelector('#position').value;
        const replacement = document.querySelector('#replacement').value;
        const mut_type = document.querySelector('#mut_type').value;
        const valued_settings = [
            '#Tm_range_min',
            '#Tm_range_max',
            '#gc_range_min',
            '#gc_range_max',
            '#length_min',
            '#length_max',
            '#flank5_range_min',
            '#flank5_range_max',
            '#flank3_range_min',
            '#flank3_range_max',
            '#forward_overlap5',
            '#forward_overlap3',
        ];
        const checked_settings = [
            '#terminate_gc',
            '#center_mutation'
        ];
        var settings = {};
        for (let set in valued_settings) {
            settings[valued_settings[set].slice(1)] = parseInt($(valued_settings[set]).val());
        }
        for (let set in checked_settings) {
            settings[checked_settings[set].slice(1)] = $(checked_settings[set]).prop('checked');
        }
        settings.primer_mode = $('#primer_mode').val();
        settings.expression_system = $('#expression_system').val();
        settings = JSON.stringify(settings);

        const loading = `
            <div class="d-flex justify-content-center">
                <div class="spinner-grow" style="width: 3rem; height: 3rem;" role="status">
                    <span class="sr-only">Loading...</span>
                </div>
            </div>
        `;
        document.querySelector('#form').innerHTML = loading;

        xhttp.open('POST', '/result');
        xhttp.onload = () => {
            const raw = xhttp.responseText;
            const data = JSON.parse(raw);
            var contents = "";
            for (let i = 0, n = Object.keys(data).length; i < n; i++) {
                contents += `
                    <div class="table-responsive text-nowrap">
                        <table class="table table-sm table-bordered table-hover">
                        <thead>
                            <tr>
                                <th scope="col"></th>
                                <th scope="col">Primer ${i + 1}</th>
                            </tr>
                        </thead>
                        <tbody>
                `;
                for (let key in data[(i+1).toString()]) {
                    contents += `
                        <tr>
                            <th scope="row">${key}</th>
                            <td>${data[(i+1).toString()][key]}</td>
                        </tr>`;
                }
                contents += `
                    </tbody></table></div>
                    `;
            }
            contents += `<button class="btn btn-blue-grey btn-rounded" onclick="window.location.href='/protein-based'">Try again</button>`;
            document.querySelector('#form').innerHTML = contents;
        };
        const data = new FormData();
        data.append('mode', 'PRO');
        data.append('sequence', sequence);
        data.append('target', target);
        data.append('position', position);
        data.append('replacement', replacement);
        data.append('mutation_type', mut_type);
        data.append('settings', settings);
        xhttp.send(data);
        return false;
    };
});
