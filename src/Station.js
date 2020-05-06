import React, { Component } from 'react';
import './App.css';
import Menu from './Menu';
import Characterize from './Characterize';
import { BrowserRouter, Route, Switch } from 'react-router-dom';
import { TransitionGroup, CSSTransition } from 'react-transition-group';
import 'mdbreact';


export default class Station extends Component {
	constructor(props) {
		super(props);
		this.state = {
			mode: null,
			pageId: 0,
			res: null,
		};

		this.changeView = this.changeView.bind(this);
		this.responseCatcher = this.responseCatcher.bind(this);
	}

	changeView(e, pageId) {
		this.setState({ pageId });
	}

	responseCatcher(e, res, mode) {
        this.setState({ res: res, mode: mode });
    }

	render() {
		return (
			<div className='container'>
				<div className='jumbotron my-5 border border-light' style={styles.appContainer}>
					<BrowserRouter>
						<Switch>
							<Route path='/characterize'>
								<Characterize responseCatcher={this.responseCatcher} changeView={this.changeView} />
							</Route>
							<Route exact path='/'>
								<Menu stations={app_data} changeView={this.changeView} />
							</Route>
						</Switch>
					</BrowserRouter>
				</div>
			</div>
		);
	}
}


const styles = {
	appContainer: {
		boxShadow: 'none',
		overflow: 'hidden',
	},
};

var app_data = [
    {
        key: 0,
        name: 'Characterization',
        publicId: 'primerdriver/characterize',
        href: '/characterize',
        color: 'danger',
    },
    {
        key: 1,
        name: 'DNA-based',
        publicId: 'primerdriver/dna',
        href: '/dna',
        color: 'primary',
    },
    {
        key: 2,
        name: 'Protein-based',
        publicId: 'primerdriver/protein.png',
        href: '/protein',
        color: 'success',
    },
];
