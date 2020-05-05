class CharacterizeView extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            cursorPosition: 1,
            sequenceLength: 0,
            submitValid: false,
            loading: false,
            mode: 'CHAR',
            sequence: '',
            mismatched_bases: 0,
            mutation_type: '',
        };
        this.formResetDefaults = { ...this.state };
    }

    textAreaHandler = (e) => {
        sequence = e.target.value
            .toUpperCase()
            .split("");
        filteredSequence = [];
        sequence.forEach((char, i) => {
            if (['A', 'T', 'C', 'G'].includes(char)) filteredSequence.push(char);
        });
        this.setState({
            sequence: filteredSequence.join(""),
            cursorPosition: e.target.selectionStart + 1,
            sequenceLength: e.target.value.length,
        });
    }

    mismatchedBasesHandler = (e) => {
        this.setState({ mismatched_bases: parseInt(e.target.value) })
    }

    mutationTypeHandler = (e) => {
        this.setState({ mutation_type: e.target.value })
    }

    resetForm = () => {
        this.setState({ ...this.formResetDefaults })
    }

    formValidator = () => {
        validSequence = (this.state.sequenceLength > 0);
        validMismatch = (this.state.mismatched_bases > 0);
        (validSequence && validMismatch) ?
            this.setState({ submitValid: true })
            :
            this.setState({ submitValid: false });
        this.formData = {
            mode: this.state.mode,
            sequence: this.state.sequence,
            mismatched_bases: this.state.mismatched_bases,
            mutation_type: this.state.mutation_type,
        }
    }

    submitHandler = (e) => {
        this.setState({ loading: true });
        xhttp = new XMLHttpRequest();
        xhttp.open('POST', '/api');
        xhttp.onload = () => {
            if (xhttp.status === 200) {
                res = JSON.parse(xhttp.responseText);
            } else {
                res = 'Request failed. Please try again later.'
            }
            this.props.responseCatcher(e, res, this.state.mode);
            this.props.changeView(e, 4);
        };
        data = new FormData();
        data.append('mode', this.state.mode);
        data.append('sequence', this.state.sequence);
        data.append('mismatched_bases', this.state.mismatched_bases);
        data.append('mutation_type', this.state.mutation_type);
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
                        className='btn btn-blue-grey btn-rounded mb-4'
                        href='/'
                        id='back'
                        onClick={(e) => this.props.changeView(e, 0)}
                    >
                        <i className='fas fa-arrow-left mr-2'></i>main menu
                    </a>
                    <h2 className='text-md-center py-2 mx-md-2 h2-responsive d-md-inline'>Characterization</h2>
                    <form
                        id='form'
                        onChange={this.formValidator}
                        onSubmit={this.submitHandler}
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
                                onChange={this.textAreaHandler}
                                onKeyUp={this.textAreaHandler}
                                onMouseUp={this.textAreaHandler}
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
                        <div className='row row-cols-1 row-cols-md-2'>
                            <div className='col form-group'>
                                <label htmlFor='mismatched_bases'>Number of mismatched bases</label>
                                <input
                                    className='form-control'
                                    min='0'
                                    type='number'
                                    id='mismatched_bases'
                                    name='mismatched_bases'
                                    value={this.state.mismatched_bases}
                                    onChange={this.mismatchedBasesHandler}
                                    required={true}
                                />
                            </div>
                            <div className='col form-group'>
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
                                    <option value='' disabled>Select mutation type</option>
                                    <option value="sub">Substitution</option>
                                    <option value="ins">Insertion</option>
                                    <option value="del">Deletion</option>
                                </select>
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
