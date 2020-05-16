import React, { Component } from 'react';
import {
    MDBIcon as Icon,
    MDBTooltip as Tooltip,
    MDBRow as Row,
    MDBCol as Col,
} from 'mdbreact';
import PropTypes from 'prop-types';


export default class SequenceInput extends Component {
    static propTypes = {
        cursorPosition: PropTypes.number.isRequired,
        handleChange: PropTypes.func.isRequired,
        sequence: PropTypes.string.isRequired,
        sequenceLength: PropTypes.number.isRequired,
    };

    render() {
        return (
            <div className='form-group'>
                <label htmlFor='sequence'>Enter sequence</label>
                <Tooltip
                    domElement
                    tag='span'
                    placement='right'
                    >
                    <span>
                        <Icon fas icon='question-circle' className='ml-2' />
                    </span>
                    <span>
                        Simply type in or paste your primer amino acid sequence. Characters will automatically be filtered to show only standard one-letter amino acid IUPAC codes.
                    </span>
                </Tooltip>
                <textarea
                    className='form-control md-textarea'
                    id='sequence'
                    name='sequence'
                    type='text'
                    value={this.props.sequence}
                    onChange={this.props.handleChange}
                    onKeyUp={this.props.handleChange}
                    onMouseUp={this.props.handleChange}
                    autoFocus
                    required
                    style={{ height: 128 }}
                    />
                <Row className='row-cols-1 row-cols-md-2'>
                    <Col className='text-left'>
                        <small className='blue-grey-text'>
                            Cursor position: {this.props.cursorPosition}
                        </small>
                    </Col>
                    <Col className='text-left text-md-right'>
                        <small className='blue-grey-text'>
                            Sequence length: {this.props.sequenceLength}
                        </small>
                    </Col>
                </Row>
            </div>
        );
    }
}
