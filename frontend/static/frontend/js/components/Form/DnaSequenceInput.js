import React, { Component } from "react";
import PropTypes from "prop-types";
import { Spacer, Row, Col, Text, Textarea, Tooltip } from "@geist-ui/react";
import * as Icon from "@geist-ui/react-icons";

export class DnaSequenceInput extends Component {
  static propTypes = {
    cursorPosition: PropTypes.number.isRequired,
    handleChange: PropTypes.func.isRequired,
    sequence: PropTypes.string.isRequired,
    sequenceLength: PropTypes.number.isRequired,
  };

  render() {
    return (
      <>
        <Row>
          <Tooltip text={"Hello"} type="dark">
            Hello
          </Tooltip>
        </Row>

        <Spacer y={0.25} />

        <Row justify="center">
          <Textarea
            placeholder="Enter sequence"
            id="sequence"
            name="sequence"
            value={this.props.sequence}
            onChange={this.props.handleChange}
            onKeyUp={this.props.handleChange}
            onMouseUp={this.props.handleChange}
            autoFocus
            required
            minHeight={128}
            width="100%"
          />
        </Row>

        <Spacer y={0.25} />

        <Row>
          <Col>
            <Text small>Cursor position: {this.props.cursorPosition}</Text>
          </Col>
          <Col>
            <Text small>Sequence length: {this.props.sequenceLength}</Text>
          </Col>
        </Row>
        {/*<label htmlFor="sequence">Enter sequence</label>*/}
        {/*<Tooltip domElement tag="span" placement="right">*/}
        {/*  <span>*/}
        {/*    <Icon fas icon="question-circle" className="ml-2" />*/}
        {/*  </span>*/}
        {/*  <span>*/}
        {/*    Simply type in or paste your primer DNA sequence. Characters will automatically be filtered to show only A,*/}
        {/*    T, C, G bases, and capitalized, as per IUPAC standards.*/}
        {/*  </span>*/}
        {/*</Tooltip>*/}
        {/*<textarea*/}
        {/*  className="form-control md-textarea"*/}
        {/*  id="sequence"*/}
        {/*  name="sequence"*/}
        {/*  type="text"*/}
        {/*  value={this.props.sequence}*/}
        {/*  onChange={this.props.handleChange}*/}
        {/*  onKeyUp={this.props.handleChange}*/}
        {/*  onMouseUp={this.props.handleChange}*/}
        {/*  autoFocus*/}
        {/*  required*/}
        {/*  style={{ height: 128 }}*/}
        {/*/>*/}
        {/*<MDBRow className="row-cols-1 row-cols-md-2">*/}
        {/*  <MDBCol className="text-left">*/}
        {/*    <small className="blue-grey-text">Cursor position: {this.props.cursorPosition}</small>*/}
        {/*  </MDBCol>*/}
        {/*  <MDBCol className="text-left text-md-right">*/}
        {/*    <small className="blue-grey-text">Sequence length: {this.props.sequenceLength}</small>*/}
        {/*  </MDBCol>*/}
        {/*</MDBRow>*/}
      </>
    );
  }
}
