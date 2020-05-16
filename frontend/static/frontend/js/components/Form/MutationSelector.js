import React, { Component } from 'react';
import {
    MDBRow as Row,
    MDBCol as Col,
} from 'mdbreact';


export default class MutationType extends Component {
    render() {
        let columns = (this.props.mutation_type === 'sub')
            ? 3
            : 2;
        if (this.props.mutation_type === '') return null;
        else return (
            <Row className={`row-cols-1 row-cols-md-${columns}`}>
                {(this.props.mutation_type !== 'ins')
                    ? <Col className='form-group'>
                        <label htmlFor='target'>Target</label>
                        <input
                            type='text'
                            className='form-control'
                            name='target'
                            id='target'
                            value={this.props.target}
                            onChange={this.props.handleChange}
                            required
                        />
                    </Col>
                    : null
                }
                <Col className='form-group'>
                    <label htmlFor='target'>Position</label>
                    <input
                        type='number'
                        min='0'
                        className='form-control'
                        name='position'
                        id='position'
                        value={this.props.position}
                        onChange={this.props.handleChange}
                        required
                    />
                </Col>
                {(this.props.mutation_type !== 'del')
                    ? <Col className='form-group'>
                        <label htmlFor='target'>Replacement</label>
                        <input
                            type='text'
                            className='form-control'
                            name='replacement'
                            id='replacement'
                            min='0'
                            value={this.props.replacement}
                            onChange={this.props.handleChange}
                            required={true}
                        />
                    </Col>
                    : null
                }
            </Row>
        );
    }
}
