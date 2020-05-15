import React, { lazy } from 'react';
import { Route, Switch } from 'react-router-dom';
import Menu from './Menu';

const Characterize = lazy(() => import('./Characterize'));
const Dna = lazy(() => import('./DnaView'));
const Protein = lazy(() => import('./ProteinView'));


export default (
    <Switch>
        <Route exact path='/' component={Menu} />
        <Route path='/characterize' component={Characterize} />
        <Route path='/dna' component={Dna} />
        <Route path='/protein' component={Protein} />
    </Switch>
);
