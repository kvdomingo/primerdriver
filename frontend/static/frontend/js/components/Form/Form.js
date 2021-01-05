import React, { Component } from "react";
import { Link } from "react-router-dom";
import { Row, Text, Button, Spacer, Grid } from "@geist-ui/react";
import * as Icon from "@geist-ui/react-icons";

export class Form extends Component {
  render() {
    return (
      <>
        <Row>
          <Link to="/" id="back">
            <Button type="secondary-light" icon={<Icon.ChevronLeftCircle />}>
              main menu
            </Button>
          </Link>
          <Spacer x={1} />
          <Text h2>{this.props.title}</Text>
        </Row>

        <Spacer y={1} />

        <Row justify="center">
          <form
            id="form"
            className="form"
            onChange={this.props.handleValidate}
            onMouseUp={this.props.handleValidate}
            onKeyUp={this.props.handleValidate}
            onSubmit={this.props.handleSubmit}
            style={{ width: "100%" }}
          >
            {this.props.children}

            <Spacer y={1} />

            <Row>
              <Button type="warning" htmlType="reset" id="reset" onClick={this.props.handleReset}>
                Reset
              </Button>
              <Spacer x={0.5} />
              <Button
                type="success"
                htmlType="submit"
                id="submit"
                onClick={this.props.handleSubmit}
                disabled={!this.props.isValid}
              >
                Submit
              </Button>
            </Row>
          </form>
        </Row>
      </>
    );
  }
}
