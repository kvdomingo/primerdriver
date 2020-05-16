import React, { Component } from 'react';
import { withRouter } from 'react-router-dom';
import LoadingScreen from './LoadingScreen';
import Form from './Form/Form';
import SequenceInput from './Form/DnaSequenceInput';
import NumberMismatch from './Form/NumberMismatch';
import MutationType from './Form/MutationType';
import {
	MDBRow as Row,
	MDBCol as Col,
} from 'mdbreact';


export default withRouter(class Characterize extends Component {
	constructor(props) {
		super(props);

		this.state = {
			formData: {},
			sequence: '',
            mismatched_bases: 0,
            mutation_type: '',
			cursorPosition: 1,
			sequenceLength: 0,
            isValid: false,
			loading: false,
            mode: 'CHAR',
        };
        this.formDefaults = { ...this.state };

		this.handleChange = this.handleChange.bind(this);
		this.handleChangeInt = this.handleChangeInt.bind(this);
		this.handleReset = this.handleReset.bind(this);
		this.handleSubmit = this.handleSubmit.bind(this);
		this.validateForm = this.validateForm.bind(this);
	}

	handleChange(e) {
		if (e.target.name === 'sequence') {
			let sequence = e.target.value
	            .toUpperCase()
	            .split("");
	        let filteredSequence = [];
	        sequence.forEach((char, i) => {
	            if (['A', 'G', 'T', 'C'].includes(char)) filteredSequence.push(char);
	        });
	        this.setState({
	            sequence: filteredSequence.join(''),
	            cursorPosition: e.target.selectionStart + 1,
	            sequenceLength: e.target.value.length,
	        });
		} else {
			let { name, value } = e.target;
			this.setState({ [name]: value });
		}
	}

	handleChangeInt(e) {
        let { name, value } = e.target;
        this.setState({ [name]: parseInt(value) })
    }

	handleReset(e) {
		this.setState({ ...this.formDefaults });
	}

	handleSubmit(e) {
		this.setState({ loading: true });
		const xhr = new XMLHttpRequest();
		xhr.open('POST', '/api');
		xhr.onload = () => {
			if (xhr.status === 200) {
				var res = JSON.parse(xhr.responseText);
			} else {
				res = 'Request failed. Please try again later.'
			}
			this.props.responseCatcher(res, this.state.mode);
			let { history } = this.props;
			history.push('/results');
		};
		const data = new FormData();
		Object.keys(this.formData).map((key, i) => (
			data.append(key, this.formData[key])
		));
		xhr.send(data);
		e.preventDefault();
	}

	validateForm(e) {
		let validSequence = (this.state.sequence.length > 0),
			validMismatch = (this.state.mismatched_bases > 0),
			validMutation = (this.state.mutation_type !== '');
		if (validSequence && validMismatch && validMutation) {
			this.setState({ isValid: true });
		} else {
			this.setState({ isValid: false });
		}
		let formData = ['mode', 'sequence', 'mismatched_bases', 'mutation_type'];
		this.formData = {};
		formData.map((item, i) => this.formData[item] = this.state[item]);
	}

	render() {
		if (this.state.loading) return (<LoadingScreen />);
		else return (
			<Form title='Primer Characterization' handleValidate={this.validateForm} handleSubmit={this.handleSubmit} handleReset={this.handleReset} {...this.state}>
				<SequenceInput handleChange={this.handleChange} {...this.state} />
				<Row className='row-cols-1 row-cols-md-2'>
					<Col>
						<NumberMismatch handleChangeInt={this.handleChangeInt} {...this.state} />
					</Col>
					<Col>
						<MutationType handleChange={this.handleChange} {...this.state} />
					</Col>
				</Row>
			</Form>
		);
	}
})
