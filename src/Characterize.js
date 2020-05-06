import React, { Component } from 'react';
import { Link } from 'react-router-dom';
import LoadingScreen from './LoadingScreen';
import 'mdbreact';


export default class Characterize extends Component {
	constructor(props) {
		super(props);
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

		this.formValidator = this.formValidator.bind(this);
		this.submitHandler = this.submitHandler.bind(this);
	}

	formValidator(e) {
		let validSequence = (this.state.sequenceLength > 0);
        let validMismatch = (this.state.mismatched_bases > 0);
        (validSequence && validMismatch) ?
            this.setState({ submitValid: true })
            :
            this.setState({ submitValid: false });
        this.formData = {
            mode: this.state.mode,
            sequence: this.state.sequence,
            mismatched_bases: this.state.mismatched_bases,
            mutation_type: this.state.mutation_type,
        };
	}

	submitHandler(e) {
		e.preventDefault();
	}

	render() {
		if (this.state.loading) return (<LoadingScreen />);
		else return (
			<div className='container mb-5'>
				<Link
					className='btn btn-blue-grey btn-rounded mb-4'
					to='/'
					id='back'
					onClick={(e) => this.props.changeView(e, 0)}
					>
					<i className='fas fa-arrow-left mr-2'></i>main menu
				</Link>
				<h2 className='text-md-center py-2 mx-md-2 h2-responsive d-md-inline'>Characterization</h2>
				<form
					id='form'
					onChange={this.formValidator}
					onSubmit={this.submitHandler}
					>
				</form>
			</div>
		);
	}
}
