import React, { Component } from "react";
import { withRouter } from "react-router-dom";
import LoadingScreen from "./LoadingScreen";
import { Form, DnaSequenceInput, NumberMismatch, MutationType } from "./Form";
import { MDBRow as Row, MDBCol as Col } from "mdbreact";

class Characterize extends Component {
  state = {
    formData: {},
    sequence: "",
    mismatched_bases: 0,
    mutation_type: "",
    cursorPosition: 1,
    sequenceLength: 0,
    isValid: false,
    loading: false,
    mode: "CHAR",
  };

  formDefaults = { ...this.state };

  handleChange = e => {
    if (e.target.name === "sequence") {
      let sequence = e.target.value.toUpperCase().split("");
      let filteredSequence = [];
      sequence.forEach(char => {
        if (["A", "G", "T", "C"].includes(char)) filteredSequence.push(char);
      });
      this.setState({
        sequence: filteredSequence.join(""),
        cursorPosition: e.target.selectionStart + 1,
        sequenceLength: e.target.value.length,
      });
    } else {
      let { name, value } = e.target;
      this.setState({ [name]: value });
    }
  };

  handleChangeInt = e => {
    let { name, value } = e.target;
    this.setState({ [name]: parseInt(value) });
  };

  handleReset = () => {
    this.setState({ ...this.formDefaults });
  };

  handleSubmit = e => {
    e.preventDefault();

    this.setState({ loading: true });

    const xhr = new XMLHttpRequest();
    xhr.open("POST", "/api");
    xhr.onload = () => {
      let res;

      if (xhr.status === 200) res = JSON.parse(xhr.responseText);
      else res = "Request failed. Please try again later.";

      this.props.responseCatcher(res, this.state.mode);
      let { history } = this.props;
      history.push("/results");
    };

    const data = new FormData();
    Object.keys(this.formData).map(key => data.append(key, this.formData[key]));
    xhr.send(data);
  };

  validateForm = () => {
    let validSequence = this.state.sequence.length > 0,
      validMismatch = this.state.mismatched_bases > 0,
      validMutation = this.state.mutation_type !== "";
    if (validSequence && validMismatch && validMutation) {
      this.setState({ isValid: true });
    } else {
      this.setState({ isValid: false });
    }
    let formData = ["mode", "sequence", "mismatched_bases", "mutation_type"];
    this.formData = {};
    formData.map(item => (this.formData[item] = this.state[item]));
  };

  render() {
    if (this.state.loading) return <LoadingScreen />;
    else
      return (
        <Form
          title="Primer Characterization"
          handleValidate={this.validateForm}
          handleSubmit={this.handleSubmit}
          handleReset={this.handleReset}
          {...this.state}
        >
          <DnaSequenceInput handleChange={this.handleChange} {...this.state} />
          <Row className="row-cols-1 row-cols-md-2">
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
}

export default withRouter(Characterize);
