import React, { lazy } from 'react';
import { Route, Switch } from 'react-router-dom';
import Menu from './Menu';
import './App.css';

const Characterize = lazy(() => import('./Characterize'));
const Dna = lazy(() => import('./DnaView'));
const Protein = lazy(() => import('./ProteinView'));
const Results = lazy(() => import('./Result'));


export default class Routes extends React.Component {
    render() {
        let { location } = this.props;
        return (
            <Switch>
                <Route exact path='/' component={Menu} />
                <Route path='/characterize' render={() => <Characterize {...this.props} />} />
                <Route path='/dna' render={() => <Dna {...this.props} />} />
                <Route path='/protein' render={() => <Protein {...this.props} />} />
                <Route path='/results' render={() => <Results results={this.props.res} mode={this.props.mode} />} />
            </Switch>
        );
    }
}
