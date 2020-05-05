class ProteinView extends React.Component {
    constructor(props) {
        super(props);
        this.state = {
            cursorPosition: 1,
            sequenceLength: 0,
            submitValid: false,
            loading: false,
            mode: 'PRO',
            sequence: '',
            mutation_type: '',
            target: '',
            position: 1,
            replacement: '',
            Tm_range_min: 75,
            Tm_range_max: 85,
            gc_range_min: 40,
            gc_range_max: 60,
            length_min: 25,
            length_max: 45,
            flank5_range_min: 11,
            flank5_range_max: 21,
            flank3_range_min: 11,
            flank3_range_max: 21,
            forward_overlap5: 9,
            forward_overlap3: 9,
            terminate_gc: true,
            center_mutation: true,
            primer_mode: 'complementary',
            expression_system: 'Homo sapiens',
        }
        this.formResetDefaults = { ...this.state }
    }

    genericNucleobaseHandler = (e) => {
        name = e.target.name;
        old_value = e.target.value
            .toUpperCase()
            .split("");
        amino_acids = ['A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        value = [];
        old_value.forEach((char, i) => {
            if (amino_acids.includes(char)) value.push(char);
        });
        value = value.join('');

        if (name === 'sequence') {
            this.setState({
                [name]: value,
                cursorPosition: e.target.selectionStart + 1,
                sequenceLength: e.target.value.length,
            });
        } else {
            this.setState({ [name]: value });
        }
    }

    mutationTypeHandler = (e) => {
        if (e.target.value !== this.state.mutation_type) {
            this.setState({
                target: '',
                replacement: '',
            });
        }
        this.setState({ mutation_type: e.target.value });
    }

    genericSelectHandler = (e) => {
        name = e.target.name;
        value = e.target.value
        this.setState({ [name]: value });
    }

    genericCheckedHandler = (e) => {
        name = e.target.name;
        value = e.target.checked
        this.setState({ [name]: value });
    }

    genericTextHandler = (e) => {
        name = e.target.name;
        value = e.target.value
        this.setState({ [name]: value });
    }

    genericFloatHandler = (e) => {
        name = e.target.name
        value = parseFloat(e.target.value)

        if ((name === 'Tm_range_min') && (value >= this.state.Tm_range_max)) return;
        if ((name === 'Tm_range_max') && (value <= this.state.Tm_range_min)) return;
        if ((name === 'gc_range_min') && (value >= this.state.gc_range_max)) return;
        if ((name === 'gc_range_max') && (value <= this.state.gc_range_min)) return;

        return this.setState({ [name]: value });
    }

    genericIntHandler = (e) => {
        name = e.target.name
        value = parseInt(e.target.value)
        this.setState({ [name]: value })
    }

    resetForm = () => {
        this.setState({ ...this.formResetDefaults });
    }

    formValidator = () => {
        validSequence = (this.state.sequenceLength > 0);
        validMutation = (this.state.mutation_type !== '');
        validMutationCode = ((this.state.target.length > 0) || (this.state.replacement.length > 0));
        (validSequence && validMutation && validMutationCode) ?
            this.setState({ submitValid: true })
            :
            this.setState({ submitValid: false });
        this.formSettings = {
            Tm_range_min: this.state.Tm_range_min,
            Tm_range_max: this.state.Tm_range_max,
            gc_range_min: this.state.gc_range_min,
            gc_range_max: this.state.gc_range_max,
            length_min: this.state.length_min,
            length_max: this.state.length_max,
            flank5_range_min: this.state.flank5_range_min,
            flank5_range_max: this.state.flank5_range_max,
            flank3_range_min: this.state.flank3_range_min,
            flank3_range_max: this.state.flank3_range_max,
            forward_overlap5: this.state.forward_overlap5,
            forward_overlap3: this.state.forward_overlap3,
            terminate_gc: this.state.terminate_gc,
            center_mutation: this.state.center_mutation,
            primer_mode: this.state.primer_mode,
            expression_system: this.state.expression_system,
        };
    }

    submitHandler = (e) => {
        this.setState({ loading: true })
        xhttp = new XMLHttpRequest();
        xhttp.open('POST', '/api');
        xhttp.onload = () => {
            if (xhttp.status === 200) {
                res = JSON.parse(xhttp.responseText);
            } else {
                res = 'Request failed. Please try again later.';
                console.log(xhttp.statusText);
                console.log(xhttp.response);
                console.log(xhttp.responseText);
            }
            this.props.responseCatcher(e, res, this.state.mode);
            this.props.changeView(e, 4);
        };
        data = new FormData();
        data.append('mode', this.state.mode);
        data.append('sequence', this.state.sequence);
        data.append('target', this.state.target);
        data.append('position', this.state.position);
        data.append('replacement', this.state.replacement);
        data.append('mutation_type', this.state.mutation_type);
        data.append('settings', JSON.stringify(this.formSettings));
        xhttp.send(data);
        e.preventDefault();
    }

    render() {
        if (this.state.loading) {
            return (
                <div className='container mb-5 text-center'>
                    <div className='spinner-grow mb-2' role='status'></div><br />
                    Please wait...
                </div>
            );
        } else {
            return (
                <div className='container mb-5'>
                    <a
                        className='btn btn-blue-grey mb-4 mr-3'
                        href='/'
                        id='back'
                        onClick={(e) => this.props.changeView(e, 0)}
                    >
                        <i className='fas fa-arrow-left mr-3'></i>main menu
                    </a>
                    <h2 className='text-md-center py-2 mx-md-2 h2-responsive d-md-inline'>Protein-based Primer Design</h2>
                    <form
                        id='form'
                        onChange={this.formValidator}
                        onKeyUp={this.formValidator}
                        onSubmit={this.submitHandler}
                        autoComplete='off'
                    >
                        <div className='form-group'>
                            <label htmlFor='sequence'>Enter DNA sequence</label>
                            <i
                                className='fas fa-question-circle ml-1'
                                data-toggle='tooltip'
                                title='Simply type in or paste your primer DNA sequence. Characters will automatically be filtered to show only A, T, C, G bases, and capitalized, as per IUPAC standards.'
                            />
                            <textarea
                                className='form-control md-textarea'
                                id='sequence'
                                name='sequence'
                                type='text'
                                value={this.state.sequence}
                                onChange={this.genericNucleobaseHandler}
                                onKeyUp={this.genericNucleobaseHandler}
                                onMouseUp={this.genericNucleobaseHandler}
                                autoFocus={true}
                                required={true}
                                style={{ height: '128px' }}
                            />
                            <div className='row row-cols-1 row-cols-md-2'>
                                <div className='col text-left'>
                                    <small className='blue-grey-text'>
                                        Cursor position: {this.state.cursorPosition}
                                    </small>
                                </div>
                                <div className='col text-md-right'>
                                    <small className='blue-grey-text'>
                                        Sequence length: {this.state.sequenceLength}
                                    </small>
                                </div>
                            </div>
                        </div>
                        <div className='form-group'>
                            <label htmlFor='mutation_type'>Mutation type</label>
                            <select
                                className='browser-default custom-select'
                                id='mutation_type'
                                name='mutation_type'
                                onChange={this.mutationTypeHandler}
                                placeholder='Select mutation type'
                                required={true}
                                value={this.state.mutation_type}
                            >
                                <option value=''>Select mutation type</option>
                                <option value="sub">Substitution</option>
                                <option value="ins">Insertion</option>
                                <option value="del">Deletion</option>
                            </select>
                        </div>
                        <MutationCodeLayout
                            { ...this.state }
                            genericNucleobaseHandler={this.genericNucleobaseHandler}
                            genericIntHandler={this.genericIntHandler}
                        />
                        <div className="accordion my-3" id="accordionAdvanced">
                            <div className="card z-depth-0 bordered">
                                <div className="card-header" id="headingAdvanced">
                                    <h5 className="mb-0 text-center">
                                        <button className="btn btn-link" id='show-advanced' type="button" data-toggle="collapse" data-target="#advanced" aria-expanded="true" aria-controls="advanced">
                                            Advanced settings
                                        </button>
                                    </h5>
                                </div>
                                <div id="advanced" className="collapse" aria-labelledby="headingAdvanced" data-parent="#accordionAdvanced">
                                    <div className="card-body">
                                        <div className="row row-cols-3 text-center">
                                            <div className="col"></div>
                                            <div className="col">Min</div>
                                            <div className="col">Max</div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                Melting point (<sup>o</sup>C)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="Tm_range_min"
                                                        name="Tm_range_min"
                                                        className="form-control"
                                                        value={this.state.Tm_range_min}
                                                        onChange={this.genericFloatHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="Tm_range_max"
                                                        name="Tm_range_max"
                                                        className="form-control"
                                                        value={this.state.Tm_range_max}
                                                        onChange={this.genericFloatHandler}
                                                    />
                                                </div>
                                            </div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                GC Content (%)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="gc_range_min"
                                                        name="gc_range_min"
                                                        className="form-control"
                                                        value={this.state.gc_range_min}
                                                        onChange={this.genericFloatHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="gc_range_max"
                                                        name="gc_range_max"
                                                        className="form-control"
                                                        value={this.state.gc_range_max}
                                                        onChange={this.genericFloatHandler}
                                                    />
                                                </div>
                                            </div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                Length (bp)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="length_min"
                                                        name="length_min"
                                                        className="form-control"
                                                        value={this.state.length_min}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="length_max"
                                                        name="length_max"
                                                        className="form-control"
                                                        value={this.state.length_max}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                5&apos; flanking region (bp)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="flank5_range_min"
                                                        name="flank5_range_min"
                                                        className="form-control"
                                                        value={this.state.flank5_range_min}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="flank5_range_max"
                                                        name="flank5_range_max"
                                                        className="form-control"
                                                        value={this.state.flank5_range_max}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                        </div>
                                        <div className="row">
                                            <div className="col text-right">
                                                3&apos; flanking region (bp)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="flank3_range_min"
                                                        name="flank3_range_min"
                                                        className="form-control"
                                                        value={this.state.flank3_range_min}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="flank3_range_max"
                                                        name="flank3_range_max"
                                                        className="form-control"
                                                        value={this.state.flank3_range_max}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                5&apos; forward overlap (bp)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="forward_overlap5"
                                                        name="forward_overlap5"
                                                        className="form-control"
                                                        value={this.state.forward_overlap5}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col"></div>
                                        </div>
                                        <div className="row row-cols-3">
                                            <div className="col text-right">
                                                3&apos; forward overlap (bp)
                                            </div>
                                            <div className="col">
                                                <div className="md-outline">
                                                    <input
                                                        type="number"
                                                        id="forward_overlap3"
                                                        name="forward_overlap3"
                                                        className="form-control"
                                                        value={this.state.forward_overlap3}
                                                        onChange={this.genericIntHandler}
                                                    />
                                                </div>
                                            </div>
                                            <div className="col"></div>
                                        </div>
                                        <div className='row row-cols-3'>
                                            <div className='col text-right'>Expression system</div>
                                            <div className='col'>
                                                <select
                                                    className='browser-default custom-select'
                                                    id='expression_system'
                                                    name='expression_system'
                                                    value={this.state.expression_system}
                                                    onChange={this.genericSelectHandler}
                                                >
                                                    {this.props.expressionList.map((exp, i) =>
                                                        <option value={exp}>{exp}</option>
                                                    )}
                                                </select>
                                            </div>
                                        </div>
                                        <div className='row row-cols-3'>
                                            <div className='col text-right'>Primer type</div>
                                            <div className='col'>
                                                <select
                                                    className="browser-default custom-select"
                                                    id="primer_mode"
                                                    name='primer_mode'
                                                    value={this.state.primer_mode}
                                                    onChange={this.genericSelectHandler}
                                                >
                                                    <option value="complementary">Complementary</option>
                                                    <option value="overlapping">Overlapping</option>
                                                </select>
                                            </div>
                                            <div className='col'></div>
                                        </div>
                                        <div className="custom-control custom-checkbox text-left mt-3">
                                            <input
                                                type="checkbox"
                                                className="custom-control-input"
                                                id="terminate_gc"
                                                name='terminate_gc'
                                                checked={this.state.terminate_gc}
                                                onChange={this.genericCheckedHandler}
                                            />
                                            <label className="custom-control-label" htmlFor="terminate_gc">Terminates in G/C</label>
                                        </div>
                                        <div className="custom-control custom-checkbox text-left">
                                            <input
                                                type="checkbox"
                                                className="custom-control-input"
                                                id="center_mutation"
                                                name='center_mutation'
                                                checked={this.state.center_mutation}
                                                onChange={this.genericCheckedHandler}
                                            />
                                            <label className="custom-control-label" htmlFor="center_mutation">Mutation at center of primer</label>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div className='text-center mt-3'>
                            <input
                                type="reset"
                                id='reset'
                                onClick={this.resetForm}
                                className="btn btn-warning text-dark"
                                value="Reset"
                            />
                            <input
                                type="submit"
                                className="btn btn-primary"
                                value="Submit"
                                id="submit"
                                disabled={!this.state.submitValid}
                            />
                        </div>
                    </form>
                </div>
            );
        }
    }
}
