import React, { Component, Suspense } from "react";
import "../App.css";
import Loading from "./LoadingScreen";
import Routes from "./Routes";
import { BrowserRouter as Router } from "react-router-dom";
import { MDBContainer as Container, MDBJumbotron as Jumbotron } from "mdbreact";

export default class Station extends Component {
  constructor(props) {
    super(props);
    this.state = {
      res: [],
      mode: "",
    };

    this.responseCatcher = this.responseCatcher.bind(this);
  }

  responseCatcher(res, mode) {
    this.setState({ res, mode });
  }

  render() {
    return (
      <Container>
        <Jumbotron className="my-5 px-md-5 border border-light" style={styles.appContainer}>
          <Suspense fallback={<Loading />}>
            <Router>
              <Routes responseCatcher={this.responseCatcher} {...this.state} />
            </Router>
          </Suspense>
        </Jumbotron>
      </Container>
    );
  }
}

const styles = {
  appContainer: {
    boxShadow: "none",
    overflow: "hidden",
  },
};
